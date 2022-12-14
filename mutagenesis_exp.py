import os
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from joblib import delayed
from joblib import Parallel


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
            mutagenesis_list.append(f"python run.py {complex_file} {amino_acid_1_letter}B{str_list[1]}{amino_acid} HL_B")

    # run job
    n_jobs = int(os.cpu_count())
    # run multiple jobs
    with Parallel(n_jobs=n_jobs, verbose=1) as parallel:
        result = parallel(
            delayed(run_geoppi)(args)
            for args in mutagenesis_list
        )

    result = np.array([np.float(x) for x in result])
    result_reshaped = result.reshape(20, len(ligand_residues))

    # create mutagenesis_exp dict:
    df = pd.DataFrame(columns=ligand_residues, data=result_reshaped, index=INTEGER_TO_RESIDUE_ONE_LETTER.tolist())
    df.style.set_precision(2).background_gradient(cmap="RdBu").hide_index().to_excel(mutagenesis_report_file, engine='openpyxl')

def ucl_project():

    # YTH24
    antibody_name = "YTH24"
    job_id = "ucl_project_part_2"
    model_name_epitope_4 = ["rank0_model1_mdref_2", "rank1_model0_mdref_23", "rank2_model1_mdref_1", "rank3_model3_mdref_20"]
    model_name_epitope_0 = ["rank12_model0_mdref_15", "rank17_model0_mdref_136", "rank22_model1_mdref_44", "rank23_model1_mdref_36"]

    for i, model in enumerate(model_name_epitope_4):
        print(f"\nRunning mutagenesis experiment for: Epitope_4: Number {i}/4: {model} ")
        mutagenesis_report(antibody_name, model, job_id)

    for i, model in enumerate(model_name_epitope_0):
        print(f"\nRunning mutagenesis experiment for: Epitope_0: Number {i}/4: {model} ")
        mutagenesis_report(antibody_name, model, job_id)

    # # YTH54
    # antibody_name = "YTH54"
    # job_id = "ucl_project_part_2"
    # model_name = ["rank4_model1_mdref_48", "rank5_model1_mdref_27"]
    #
    # for model in model_name:
    #     print("\nRunning mutagenesis experiment for:", model)
    #     mutagenesis_report(antibody_name, model, job_id)

ucl_project()