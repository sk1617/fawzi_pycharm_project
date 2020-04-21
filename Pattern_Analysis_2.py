import glob
import pandas as pd
import numpy as np
import statistics as st
import matplotlib.pyplot as plt
import re
from DataProcessing import *
from DataImporting import sequence_list
from venditti_data_importer import venditti_assignment
peak_list, residue_list = main_data_processing()

acceptance_threshold = .20
should_print_unassigned_lists = False
should_print_peak_data = True
should_export_to_excel = True
should_take_best_10 = False


# if you want to analyze multiple test conditions, add them here.
for random_index in range(1):
    # DATA INGESTION
    time_taken_list = list()
    final_energy_list = list()
    filename_list = list()
    count_of_assignments = 0

    # glob.glob returns a list of paths that match filename NAME.
    # add all files for each test condition

    # NAME = ["Slurm Trials/slurm-4102660_{}.out".format(str(x)) for x in range(1, 17)] + \
    #        ["Slurm Trials/slurm-4109433_{}.out".format(str(x)) for x in range(17, 201)] + \
    #        ["Slurm Trials/slurm-9199622_{}.out".format(str(x)) for x in range(1, 81)]

    NAME = ["mmals/Fawzi_pycharm_project/slurm-new_and_improved?_{}.out".format(str(x)) for x in range(1, 200)]



    set_iterations, initial_temperature, exponent, eif = None, None, None, None
    time_taken, completed_iterations, final_temp, og_index_list, og_energy, index_list, final_energy = \
        None, None, None, None, None, None, None

    initial_data = pd.DataFrame(columns=['set_iterations', 'time_taken', 'completed_iterations', 'final_temp',
                                         'og_index_list', 'og_energy', 'index_list', 'final_energy'])

    for filename_str in NAME:
        # imports each file
        # filename = [i[0] for i in glob.glob(filename_str)]
        try:
            file = open(filename_str, 'r')
            filename_list.append(filename_str)


            for line in file.readlines():
                if "Done" in line:
                    count_of_assignments += 1
                    df_temp = pd.DataFrame([[set_iterations, time_taken, completed_iterations,
                                             final_temp, og_index_list, og_energy,
                                             index_list, final_energy]],
                                           columns=['set_iterations', 'time_taken', 'completed_iterations',
                                                    'final_temp', 'og_index_list', 'og_energy',
                                                    'index_list', 'final_energy'])
                    initial_data = initial_data.append(df_temp, ignore_index=True)

                    set_iterations, initial_temperature, exponent, eif = None, None, None, None
                    time_taken, completed_iterations, final_temp, og_index_list, og_energy, index_list, final_energy = \
                        None, None, None, None, None, None, None

                elif 'exponent:' in line:
                    set_iterations = int(re.search(r" \d* ", line)[0][1:-1])
                elif 'time taken (sec)' in line:
                    time_taken = float(re.search(r"\d*\.\d*", line)[0])
                    time_taken_list.append(time_taken)
                elif 'number of iterations' in line:
                    completed_iterations = int(re.search(r"[0-9].*", line)[0])
                elif 'final temperature' in line:
                    final_temp = float(re.search(r"\d*\.\d*", line)[0])
                    final_energy_list.append(final_temp)
                elif 'og index list:' in line:
                    og_index_list = list()
                    [og_index_list.append(int(x)) for x in line[line.index('[') + 1:line.index(']')].split(', ')]
                elif 'og index list energy:' in line:
                    og_energy = float(re.search(r"\d*\.\d*", line)[0])
                elif 'index list:' in line:
                    index_list = list()
                    [index_list.append(int(x)) for x in line[line.index('[') + 1:line.index(']')].split(', ')]
                elif 'index list energy:' in line:
                    final_energy = float(re.search(r"\d*\.\d*", line)[0])
        except FileNotFoundError:
            print(filename_str)

    initial_data = initial_data.sort_values(by='final_energy')
    print(initial_data.head(), initial_data.columns)
    data = pd.DataFrame(columns=['residue number', 'amino acid', 'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime',
                             'CACAPrime', 'CBCBPrime', 'CACBExp', 'CACBExpPrime', 'NOESY', 'NOESYSubWeights', 'NOESYDelH'])
    if True:
        index_list = initial_data.loc[158, 'index_list']
        info_list = eval_energy(index_list, peak_list, residue_list, energy_if_false)[1]
        for i, pai in enumerate(index_list):
            residue_number = i + 1
            amino_acid = sequence_list[i]

            peak = peak_list[pai]
            assert pai == peak.get_data('peakNumber')
            n_shift = peak.get_data('TROSYNShift')
            h_shift = peak.get_data('TROSYHShift')
            ca_shift = peak.get_data('CAShift')
            ca_prime_shift = peak.get_data('CAPrimeShift')
            cb_shift = peak.get_data('CBShift')
            cb_prime_shift = peak.get_data('CBPrimeShift')

            ca_ca_prime, cb_cb_prime, ca_cb_exp, = energy_if_false, energy_if_false, energy_if_false
            ca_cb_exp_prime, noesy_energy = energy_if_false, energy_if_false
            noesy_subweights, noesy_del_h = [], []

            for k, line in enumerate(info_list):
                if len(line) == 2:
                    if line[0] == i:
                        assert line[1] == pai
                        for l in range(1,6):
                            if len(info_list[k + l]) == 2:
                                break
                            elif 'CACAPrime' == info_list[k + l][0]:
                                ca_ca_prime = info_list[k + l][1]
                            elif 'CBCBPrime' == info_list[k + l][0]:
                                cb_cb_prime = info_list[k + l][1]
                            elif 'CACBExp' == info_list[k + l][0]:
                                ca_cb_exp = info_list[k + l][1]
                            elif 'CACBExpPrime' == info_list[k + l][0]:
                                ca_cb_exp_prime = info_list[k + l][1]
                            elif 'NOESY' == info_list[k + l][0]:
                                noesy_energy = info_list[k + l][1]
                                noesy_subweights = info_list[k + l][2]
                                noesy_del_h = info_list[k + l][3]

            temp_df = pd.DataFrame([[residue_number, amino_acid, n_shift, h_shift, ca_shift, ca_prime_shift, cb_shift, cb_prime_shift,
                                     ca_ca_prime, cb_cb_prime, ca_cb_exp, ca_cb_exp_prime, noesy_energy, noesy_subweights, noesy_del_h]],
                                   columns=['residue number', 'amino acid', 'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime',
                                            'CACAPrime', 'CBCBPrime', 'CACBExp', 'CACBExpPrime', 'NOESY', 'NOESYSubWeights', 'NOESYDelH'])
            data = data.append(temp_df, ignore_index=True)
        print(data.head())
        data.to_csv('/Volumes/Transcend/fawzi_pycharm_project/data_TEST_{}.csv'.format('sources_of_energy'))
