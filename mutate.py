#!/usr/bin/env python3
# -*-coding:utf-8 -*-

import os
from pathlib import Path
from subprocess import Popen, PIPE
# import openmm.app as app
# import pdbfixer
import numpy as np
import re


# paths to files
BIN_PATH = os.getcwd()
EVO_EF1_BIN = os.path.join(BIN_PATH, "EvoEF1", "EvoEF")
EVO_EF2_BIN = os.path.join(BIN_PATH, "EvoEF2", "EvoEF")


def mutate_EvoEF1(wildtype_file_name, mutant_dir, mutation_info_list):

    # write mutation information to text file (individual_list_) for input to foldx
    wildtype_file_path = Path(wildtype_file_name)
    wildtype_name = wildtype_file_path.stem
    mutant_name = wildtype_name + "_" + "_".join(mutation_info_list)

    # create unique dir name
    # when running parrallel processing, the output file "Model_001.pdb" must be unique
    os.system(f"cp {wildtype_file_name} {os.path.join(wildtype_file_path.parent, mutant_name + '.pdb')}")
    wildtype_file_path = Path(os.path.join(wildtype_file_path.parent, mutant_name + '.pdb'))


    # create mutation file input for EvoEF2
    individual_list_name = "individual_list_" + mutant_name + ".txt"
    individual_list_file = os.path.join(wildtype_file_path.parent, individual_list_name)
    with open(individual_list_file, 'w') as f:
        mutation_str_list = ""
        for mutation in mutation_info_list:
            mutation_str = '{},'.format(mutation)
            mutation_str_list += mutation_str
        mutation_str_list = re.sub(r".$", ";", mutation_str_list)
        f.write(mutation_str_list)

    # run EvoEF
    args = [EVO_EF1_BIN, "--command=BuildMutant", "--pdb",
            wildtype_file_path.name, "--mutant_file", individual_list_name,
            "[--num_of_runs=10]"]
    process = Popen(args, stdout=PIPE, stderr=PIPE, cwd=wildtype_file_path.parent)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('UTF-8')
    stderr = stderr.decode('UTF-8')
    print(stdout, stderr)

    # remove other generated file
    os.remove(os.path.join(wildtype_file_path.parent, individual_list_name))

    # change the name of the mutant file
    os.system('mv {}/{}_Model_0001.pdb   {}/{}.pdb '.format(wildtype_file_path.parent, wildtype_file_path.stem, mutant_dir,
                                                            mutant_name))

    # remove identity mutation WT
    os.system(f"rm {wildtype_file_path.parent}/{wildtype_file_path.stem}_Model_0001_WT.pdb")


def mutate_EvoEF2(wildtype_file_name, mutant_dir, mutation_info_list):
    # write mutation information to text file (individual_list_) for input to foldx
    wildtype_file_path = Path(wildtype_file_name)
    wildtype_name = wildtype_file_path.stem
    mutant_name = wildtype_name + "_" + "_".join(mutation_info_list)

    # create mutation file input for EvoEF2
    individual_list_name = "individual_list_" + mutant_name + ".txt"
    individual_list_file = os.path.join(wildtype_file_path.parent, individual_list_name)
    with open(individual_list_file, 'w') as f:
        mutation_str_list = ""
        for mutation in mutation_info_list:
            mutation_str = '{},'.format(mutation)
            mutation_str_list += mutation_str
        mutation_str_list = re.sub(r".$", ";", mutation_str_list)
        f.write(mutation_str_list)

    # run EvoEF2
    args = [EVO_EF2_BIN, "--command=BuildMutant", "--pdb",
            wildtype_file_path.name, "--mutant_file", individual_list_name]
    process = Popen(args, stdout=PIPE, stderr=PIPE, cwd=wildtype_file_path.parent)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('UTF-8')
    stderr = stderr.decode('UTF-8')
    print(stdout, stderr)

    # remove other generated file
    os.remove(os.path.join(wildtype_file_path.parent, individual_list_name))

    # change the name of the mutant file
    os.system('mv {}/{}_Model_0001.pdb   {}/{}.pdb '.format(wildtype_file_path.parent, wildtype_name, mutant_dir,
                                                            mutant_name))

