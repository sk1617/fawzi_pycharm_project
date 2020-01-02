import glob
import pandas as pd
import numpy as np
import statistics as st
import matplotlib.pyplot as plt
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
    uninitialized_table = list()
    time_taken_list = list()
    final_energy_list = list()
    filename_list = list()
    count_of_assignments = 0
    # extracts index list
    # glob.glob returns a list of paths that match filename NAME.
    # add all files for each test condition

    # NAME = ["Slurm Trials/slurm-4102660_{}.out".format(str(x)) for x in range(1, 17)] + \
    #        ["Slurm Trials/slurm-4109433_{}.out".format(str(x)) for x in range(17, 201)]

    # NAME = ["Slurm Trials/slurm-9199622_{}.out".format(str(x)) for x in range(1,81)]

    NAME = ["Slurm Trials/NOESY_80.out"]

    for filename_str in NAME:
        # imports each file
        # filename = [i[0] for i in glob.glob(filename_str)]
        file = open(filename_str, 'r')
        filename_list.append(filename_str)
        for line in file:
            pass


    assert len(uninitialized_table) > 0

    # Makes numpy array, table, and makes pandas dataframe, data
    table = np.array(uninitialized_table)
    table_inverted = table.T


    # DATA PROCESSING
    # mode table
    mode_list_full = list()
    mode_list = list()
    unassigned_positions = list()



    # Sorting data into mode list and into data dataframe
    for i, assignment_list in enumerate(table_inverted):
        assignment_list = list(assignment_list)

        # PREVIOUS
        try:
            mode = int(st.mode(assignment_list))
            count = assignment_list.count(mode)
        except st.StatisticsError:
            count = 0
        # only adds to mode list if more than threshold % of assignments agree
        if count/count_of_assignments > acceptance_threshold:
            mode_list_full.append((mode, count))
            mode_list.append(mode)
        # advanced case for Prolines and unassigned
        else:
            temporary_list = list()
            for assignment in assignment_list:
                peak = peak_list[assignment]
                peak_number = peak.get_data('peakNumber')
                trosy_h_shift = peak.get_data('TROSYHShift')
                assert assignment == peak_number
                temporary_list.append(trosy_h_shift)
            # must add exception if the number of total Null peaks > threshold
            if temporary_list.count(None)/count_of_assignments > acceptance_threshold:
                mode_list_full.append((None, temporary_list.count(None)))
                mode_list.append(None)
            else:
                mode_list_full.append(assignment_list)
                mode_list.append('NA')

    # which peaks appear twice in mode_list
    duplicated_assignments = []
    for i in range(0, len(sequence_list)):
        if mode_list.count(i) > 1:
            duplicated_assignments.append(i)


    # NEW (Adding to "data" dataframe)
    data = pd.DataFrame(columns=['residue number', 'amino acid', 'mode', 'freq',
                                 'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime', 'duplicated?', 'venditti assignment'])

    # 1. Add rows to data DataFrame
    # proline_index_list = [i for i, ltr in enumerate(sequence_list) if ltr == 'P']
    for i, line in enumerate(mode_list_full):
        residue_number = i + 1
        amino_acid = sequence_list[i]

        # find mode and freq
        if amino_acid == 'P':
            mode = None
            freq = 0
        elif type(line) is tuple:
            mode = line[0]
            freq = round(line[1]/count_of_assignments, 4)
        elif type(line) is list:
            try:
                mode = st.mode(line)
                freq = round(line.count(mode)/count_of_assignments, 4)
            except st.StatisticsError:
                mode = 'could not find mode'
                freq = 'NA'

        # find other data
        if type(mode) == int:
            peak = peak_list[mode]
            assert mode == peak.get_data('peakNumber')
            n_shift = peak.get_data('TROSYNShift')
            h_shift = peak.get_data('TROSYHShift')
            ca_shift = peak.get_data('CAShift')
            ca_prime_shift = peak.get_data('CAPrimeShift')
            cb_shift = peak.get_data('CBShift')
            cb_prime_shift = peak.get_data('CBPrimeShift')
            duplicated = True if mode in duplicated_assignments else False
        else:
            n_shift = None
            h_shift = None
            ca_shift = None
            ca_prime_shift = None
            cb_shift = None
            cb_prime_shift = None
            duplicated = False

        venditti = venditti_assignment[i]


        temp_df = pd.DataFrame([[residue_number, amino_acid, mode, freq,
                                 n_shift, h_shift, ca_shift, ca_prime_shift, cb_shift, cb_prime_shift, duplicated, venditti]],
                               columns=['residue number', 'amino acid', 'mode', 'freq',
                                        'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime', 'duplicated?', 'venditti assignment'])
        data = data.append(temp_df, ignore_index=True)

    # Data visualization

    # DATA ANALYSIS
    sequence_list = sequence_list
    len_index_list = len(sequence_list)

    # figure out how many non-Null peaks there should be
    count_non_null_peaks = len(sequence_list) - sequence_list.count('P')

    # which peaks do not appear
    not_assigned = []
    for i in range(0, count_non_null_peaks):
        if i not in mode_list:
            not_assigned.append(i)

    print('the following files were analyzed: {} with a total of {} assignments'.format(filename_list, count_of_assignments))
    print('the average time taken for each trial was: {} sec'.format(st.mean(time_taken_list)))
    print('the average end energy was: {}'.format(st.mean(final_energy_list)))
    print('duplicate peaks: {}'.format(duplicated_assignments))
    print('not assigned peaks: {}'.format(not_assigned))

    for i, line in enumerate(mode_list_full):
        # if successful assignment
        if type(line) == tuple:
            print(sequence_list[i], ':', 'mode:', line[0],
                          'count:{} ({}%)'.format(line[1], round(100 * line[1]/count_of_assignments, 5)))

        # if not successful and not None
        elif type(line) == list:
            try:
                if should_print_unassigned_lists:
                    print(sequence_list[i], ':', line, 'mode:', st.mode(line), 'count:', line.count(st.mode(line)))
                else:
                    print(sequence_list[i], ':', 'mode:', st.mode(line),
                          'count:{} ({}%)'.format(line.count(st.mode(line)),
                                                  round(100 * line.count(st.mode(line))/count_of_assignments, 5)))
            except st.StatisticsError:
                if should_print_unassigned_lists:
                    print(sequence_list[i], ':', line.sort())
                else:
                    print(sequence_list[i], ':', "No mode")
        # if none
        else:
            print(sequence_list[i], ':', line)
    print('DONE\n\n')

    if should_export_to_excel: data.to_csv('/Volumes/Transcend/Fawzi_pycharm_project/data_long.csv')
