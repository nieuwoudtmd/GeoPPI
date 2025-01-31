import json
import os
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from joblib import delayed
from joblib import Parallel
import time
from mutate import mutate_EvoEF1

INTEGER_TO_RESIDUE_ONE_LETTER = np.array([
    "A", "C", "D", "E", "F",
    "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y"])

INTEGER_TO_RESIDUE_THREE_LETTER = [
    "ALA", "CYS", "ASP", "GLU", "PHE",
    "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG",
    "SER", "THR", "VAL", "TRP", "TYR"]


def run_geoppi(command):
    # run arpeggio
    args = command.split(" ")
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('UTF-8')
    stderr = stderr.decode('UTF-8')

    return stdout


class Mutation:
    def __init__(self, pdb_dir: str, mutation_info_list: list):
        """
        pdb_dir (str): path to pdb directory with mutant and wildtype file
        mutation_info_list (list): list containing mutation information
        """

        # assign list information to variables
        self.pdb_id = mutation_info_list[0]
        self.mutation_chain = mutation_info_list[1]
        self.mutation_residue_index = mutation_info_list[2]
        self.receptor_chain = mutation_info_list[3]
        self.ligand_chain = mutation_info_list[4]
        self.mutation_info = mutation_info_list[5]

        # path to files
        self.wildtype_file = os.path.join(pdb_dir, f"{self.pdb_id}_{self.receptor_chain}_{self.ligand_chain}.pdb")
        self.mutant_file = os.path.join(pdb_dir,
                                        f"{self.pdb_id}_{self.receptor_chain}_{self.ligand_chain}_{self.mutation_info}.pdb")

        # check if file exist
        if not os.path.isfile(self.wildtype_file):
            print(f"ERROR: CANNOT FIND WILDTYPE PDB: {self.wildtype_file}")


class Mutagenesis:
    def __init__(self, project_dir, job_id):
        self.project_dir = project_dir
        self.job_id = job_id

        # directory paths
        self.job_dir = os.path.join(project_dir, job_id)
        self.pdb_dir = os.path.join(self.job_dir, "pdb")
        self.embedding_dir = os.path.join(self.job_dir, "embeddings")
        self.list_dir = os.path.join(self.job_dir, "lists")
        os.makedirs(self.pdb_dir, exist_ok=True)
        os.makedirs(self.embedding_dir, exist_ok=True)
        os.makedirs(self.list_dir, exist_ok=True)

        # file path
        self.mutation_entries_file = os.path.join(self.list_dir, "entries.json")

        # read mutation info file
        self.mutation_entries = json.load(open(self.mutation_entries_file))

    def generate_mutants(self):

        complex_file_list = []
        mutant_dir_list = []
        mutation_info_list = []
        for mutation_item in self.mutation_entries:

            mutation = Mutation(self.pdb_dir, mutation_item)

            if not os.path.isfile(os.path.join(self.pdb_dir, mutation.mutant_file)):
                complex_file_list.append(mutation.wildtype_file)
                mutant_dir_list.append(self.pdb_dir)
                mutation_info_list.append([mutation.mutation_info])

        # run job
        n_jobs = int(os.cpu_count())
        # run multiple jobs
        with Parallel(n_jobs=n_jobs, verbose=1) as parallel:
            parallel(
                delayed(mutate_EvoEF1)(complex_file, mutant_dir, mutation_info)
                for complex_file, mutant_dir, mutation_info in
                zip(complex_file_list, mutant_dir_list, mutation_info_list)
            )

    def run_geoppi(self):
        pass

    def run_geoppi_parallel(self, run_script):

        command_line_list = []
        for mutation_item in self.mutation_entries:
            # create Mutation object
            mutation = Mutation(self.pdb_dir, mutation_item)
            command_line_list.append(
                f"python {run_script}.py {mutation.wildtype_file} {mutation.mutation_info} {mutation.receptor_chain}_{mutation.ligand_chain} {self.job_dir}")

        # run job
        n_jobs = int(os.cpu_count())
        # run multiple jobs
        with Parallel(n_jobs=n_jobs, verbose=1) as parallel:
            result = parallel(
                delayed(run_geoppi)(args)
                for args in command_line_list
            )
        geoppi_result = np.array([float(x) for x in result])

        # save in mutation_list_file
        self.save_results(geoppi_result)

    def save_results(self, results):

        for mutation_item, prediction in zip(self.mutation_entries, results):
            mutation_item[6] = prediction

        mutation_list_geoppi = os.path.join(self.list_dir, "entries_geoppi_predictions.json")
        with open(mutation_list_geoppi, 'w') as f:
            json.dump(self.mutation_entries, f)


def mutagenesis_report(antibody_name, model_name, job_id):
    # paths
    data_dir = os.path.join("data", job_id, antibody_name)
    output_dir = os.path.join("mutagenesis_jobs", job_id, antibody_name)
    os.makedirs(output_dir, exist_ok=True)

    # path to pdb file
    complex_file = os.path.join(data_dir, "pdb", model_name + ".pdb")
    molecular_interaction_report_file = os.path.join(data_dir, "molecular_interaction_report", model_name + ".xlsx")
    mutagenesis_report_file = os.path.join(output_dir, model_name + ".xlsx")

    if os.path.isfile(mutagenesis_report_file):
        return

    # get interface residues
    # read molecular interaction
    df = pd.read_excel(molecular_interaction_report_file, sheet_name="residue_interaction")

    # get unique list of
    ligand_residues = df["Ligand Residue"].unique()

    mutagenesis_list = []
    for ligand_residue in ligand_residues:
        for amino_acid in INTEGER_TO_RESIDUE_ONE_LETTER:
            str_list = ligand_residue.split("/")
            amino_acid_idx = np.where(np.array(INTEGER_TO_RESIDUE_THREE_LETTER) == np.array(str_list[2]))
            amino_acid_1_letter = INTEGER_TO_RESIDUE_ONE_LETTER[amino_acid_idx[0]][0]
            mutagenesis_list.append(
                f"python run.py {complex_file} {amino_acid_1_letter}B{str_list[1]}{amino_acid} HL_B")

    # run job
    n_jobs = int(os.cpu_count())
    # run multiple jobs
    with Parallel(n_jobs=n_jobs, verbose=1) as parallel:
        result = parallel(
            delayed(run_geoppi)(args)
            for args in mutagenesis_list
        )

    result = np.array([np.float(x) for x in result])
    result_reshaped = result.reshape(len(ligand_residues), 20).T

    # # result test
    # results = np.random.uniform(low=-2, high=2, size=(20*len(ligand_residues),))
    # resutls_reshaped = results.reshape(len(ligand_residues), 20)
    # resutls_reshaped_transposed = resutls_reshaped.T

    # create mutagenesis_exp dict:
    df = pd.DataFrame(columns=ligand_residues, data=result_reshaped, index=INTEGER_TO_RESIDUE_ONE_LETTER.tolist())
    df.style.set_precision(2).background_gradient(cmap="RdBu").hide_index().to_excel(mutagenesis_report_file,
                                                                                     engine='openpyxl')


def ucl_project():
    project_dir = os.path.join(os.getcwd(), "data", "ucl_CD40")

    # antibody_name = "YTH54"
    antibody_name = "YTH24"

    # antibody related epitope dict
    epitope_hotspot_dict = {"YTH24": {
        "epitope_4": ["rank0_model1_mdref_2", "rank1_model0_mdref_23", "rank2_model1_mdref_1", "rank3_model3_mdref_20"],
        "epitope_0": ["rank12_model0_mdref_15", "rank17_model0_mdref_136", "rank22_model1_mdref_44",
                      "rank23_model1_mdref_36"]},
        "YTH54": {
            "epitope_1_1": ["rank4_model1_mdref_48", "rank30_model2_mdref_60", "rank36_model2_mdref_76"],
            "epitope_1_2": ["rank5_model1_mdref_27", "rank11_model2_mdref_53", "rank18_model2_mdref_54",
                            "rank24_model1_mdref_29"]}}

    # epitope hotspots for antibody
    epitope_complex_list = epitope_hotspot_dict[antibody_name]["epitope_0"]

    for i, experiment in enumerate(epitope_complex_list):
        time_start = time.time()
        print(f"Running job({i + 1}/{len(epitope_complex_list)}): {experiment} ")
        mutagenesis_exp = Mutagenesis(project_dir, experiment)
        mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
        print(f"Experiment run time: {time_start - time.time()}")


def bayer_project(project_dir, epitope_complex_list):
    for i, experiment in enumerate(epitope_complex_list):
        time_start = time.time()
        print(f"Running job({i + 1}/{len(epitope_complex_list)}): {experiment} ")
        mutagenesis_exp = Mutagenesis(project_dir, experiment)
        mutagenesis_exp.generate_mutants()
        mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
        print(f"Experiment run time: {time_start - time.time()}")


def mutagenesis_experiment():
    # complex list to evaluate in study
    complex_list = ["3HFM", "3NGB", "1MHP", "1JRH", "1VFB", "2JEL"]
    complex_list = ["1DQJ", "1DVF", "1CHO", "1PPF"]
    complex_list = ["3SGB", "1R0R", "3BN9", "1CBW", "3C60"]
    complex_list = ["3HFM", "3NGB", "1MHP", "1JRH"]

    # epitope setup
    epitope_list = ["epitope_original", "epitope_top_1", "epitope_top_2", "epitope_diff_1", "epitope_diff_2",
                    "epitope_diff_3"]
    epitope_list = ["epitope_original", "epitope_top_1", "epitope_top_2", "epitope_diff_1", "epitope_diff_2",
                    "epitope_diff_3", "epitope_diff_4", "epitope_diff_5"]
    epitope_list = ["mdref_1", "mdref_2", "mdref_3"]

    for complex in complex_list:
        project_dir = os.path.join(os.getcwd(), "data", "mutagenesis_experiments", complex)
        project_dir = os.path.join(os.getcwd(), "data", "mutagenesis_experiments_ab_modelled", complex)

        for i, experiment in enumerate(epitope_list):
            time_start = time.time()
            print(f"Running job({i + 1}/{len(epitope_list)}): {experiment} ")
            mutagenesis_exp = Mutagenesis(project_dir, experiment)
            print("Generating mutants")
            mutagenesis_exp.generate_mutants()
            mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
            print(f"Experiment run time: {time_start - time.time()}")


def epitope_hunt():
    # complex list to evaluate in study
    complex_list = ["3HFM", "1MHP", "1JRH"]

    # epitope setup
    epitope_list = [f"epitope_cluster_{x}_model_{y}" for x in range(1, 21) for y in range(1, 3)]

    for complex in complex_list:
        project_dir = os.path.join(os.getcwd(), "data", "epitope_hunt", complex)

        for i, experiment in enumerate(epitope_list):
            time_start = time.time()
            print(f"Running job({i + 1}/{len(epitope_list)}): {experiment} ")
            mutagenesis_exp = Mutagenesis(project_dir, experiment)
            print("Generating mutants")
            mutagenesis_exp.generate_mutants()
            mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
            print(f"Experiment run time: {time_start - time.time()}")


def bayer_part2():
    # complex list to evaluate in studies
    experiment_list = ["rank0_model3_mdref_258_r_b_rank0_model3_mdref_258_l_b_paratope_EveParatope_epitope_no_info",
                       "rank3_model3_mdref_194_r_b_rank3_model3_mdref_194_l_b_paratope_EveParatope_epitope_no_info",
                       "rank0_model3_mdref_258_r_b_rank0_model3_mdref_258_l_b_paratope_EveParatopeANDexp_epitope_no_info",
                       "rank3_model3_mdref_194_r_b_rank3_model3_mdref_194_l_b_paratope_EveParatopeANDexp_epitope_no_info"]

    # epitope setup
    epitope_list = [f"mdref_{i}" for i in range(1, 11)]

    for experiment in experiment_list:
        project_dir = os.path.join(os.getcwd(), "data", "bayer_IL11", experiment)

        for i, experiment in enumerate(epitope_list):
            time_start = time.time()
            print(f"Running job({i + 1}/{len(epitope_list)}): {experiment} ")
            mutagenesis_exp = Mutagenesis(project_dir, experiment)
            print("Generating mutants")
            mutagenesis_exp.generate_mutants()
            mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
            print(f"Experiment run time: {time_start - time.time()}")


def mutagenesis_experiment_covid():
    # complex list to evaluate in study
    complex_list = ["C002", "C104", "C105", "C119", "C144"]

    for i, complex in enumerate(complex_list):
        project_dir = os.path.join(os.getcwd(), "data", "mutagenesis_experiment_covid")

        time_start = time.time()
        print(f"Running job({i + 1}/{len(complex_list)}): {complex} ")
        mutagenesis_exp = Mutagenesis(project_dir, complex)
        mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
        print(f"Experiment run time: {time_start - time.time()}")


def mutagenesis_experiment_covid_2():
    # complex list to evaluate in study
    complex_list = ["6M0J"]

    for i, complex in enumerate(complex_list):
        project_dir = os.path.join(os.getcwd(), "data", "mutagenesis_experiment_covid_2")

        time_start = time.time()
        print(f"Running job({i + 1}/{len(complex_list)}): {complex} ")
        mutagenesis_exp = Mutagenesis(project_dir, complex)
        mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
        print(f"Experiment run time: {time_start - time.time()}")


def bayer_PAD4():

    project_name = "bayer_PAD4"
    experiment_name = "md_simulation_bayer"
    job_id_list = [f"B276_cluster_{x}" for x in range(0, 10)]

    project_dir = os.path.join(os.getcwd(), "data", project_name, experiment_name, "experiments")

    for i, job_id in enumerate(job_id_list):
        time_start = time.time()
        print(f"Running job({i + 1}/{len(job_id_list)}): {job_id} ")
        mutagenesis_exp = Mutagenesis(project_dir, job_id)
        print("Generating mutants")
        mutagenesis_exp.generate_mutants()
        mutagenesis_exp.run_geoppi_parallel(run_script="run_silicogenesis")
        print(f"Experiment run time: {time_start - time.time()}")


bayer_PAD4()
