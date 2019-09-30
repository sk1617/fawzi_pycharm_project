import glob
import pandas as pd
import numpy as np
import statistics as st
from DataProcessing import *
from full_DataImporting import sequence_list
peak_list, residue_list = main_data_processing()

acceptance_threshold = 0.5
should_print_unassigned_lists = False

for NAME in ['MMAl.out', 'exponential.out', 'linear.out']:
    # DATA INGESTION
    uninitialized_table = list()
    time_taken_list = list()
    final_energy_list = list()
    filename_list = list()
    count_of_assignments = 0
    # extracts index list
    # gets list of file names that match '*.out'
    for filename in glob.glob(NAME):
        # imports each file
        file = open(filename, 'r')
        filename_list.append(filename)
        for line in file:
            # selects the correct line
            if 'index list' in line and not 'og index list' in line and not 'energy' in line:
                # finds [ and ] in str to made
                start_index = line.index('[') + 1
                end_index = line.index(']')
                line_list_str = line[start_index:end_index]
                # makes IL with ints
                index_list = line_list_str.split(',')
                index_list = [int(x) for x in index_list]
                # asserts there are no repeats in IL
                assert len(index_list) == len(set(index_list))
                uninitialized_table.append(index_list)
                count_of_assignments += 1
            if 'time taken (sec):' in line:
                start_index = line.index(':') + 2
                time_taken_str = line[start_index:]
                time_taken = float(time_taken_str)
                time_taken_list.append(time_taken)
            if 'index list energy:' in line and not 'og' in line:
                start_index = line.index(':') + 2
                final_energy_str = line[start_index:]
                final_energy = float(final_energy_str)
                final_energy_list.append(final_energy)

    # Makes numpy array, table, and makes pandas dataframe, data
    table = np.array(uninitialized_table)
    table_inverted = table.T
    data = pd.DataFrame(table)

    # DATA PROCESSING
    # mode table
    mode_list_full = list()
    mode_list = list()
    unassigned_positions = list()

    for i, assignment_list in enumerate(table_inverted):
        assignment_list = list(assignment_list)
        try:
            mode = int(st.mode(assignment_list))
            count = assignment_list.count(mode)
        except st.StatisticsError:
            count = 0
        # only adds to mode list if more than threshold % of assignments agree
        if count/count_of_assignments > acceptance_threshold:
            mode_list_full.append(mode)
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
                mode_list_full.append(None)
                mode_list.append(None)
            else:
                mode_list_full.append(assignment_list)
                mode_list.append('NA')

    # DATA ANALYSIS
    sequence_list = sequence_list
    len_index_list = len(sequence_list)

    # figure out how many non-Null peaks there should be
    count_non_null_peaks = len(sequence_list) - sequence_list.count('P')

    # which peaks appear twice in mode_list
    duplicated_assignments = []
    for i in range(0, len_index_list):
        if mode_list.count(i) > 1:
            duplicated_assignments.append(i)

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

    # for i, line in enumerate(mode_list_full):
    #     # if successful assignment
    #     if type(line) == int:
    #         print(sequence_list[i], ':', line)
    #     # if not successful and not None
    #     elif type(line) == list:
    #         try:
    #             if should_print_unassigned_lists:
    #                 print(sequence_list[i], ':', line, 'mode:', st.mode(line), 'count:', line.count(st.mode(line)))
    #             else:
    #                 print(sequence_list[i], ':', 'mode:', st.mode(line), 'count:', line.count(st.mode(line)))
    #         except st.StatisticsError:
    #             if should_print_unassigned_lists:
    #                 print(sequence_list[i], ':', line.sort())
    #             else:
    #                 print(sequence_list[i], ':', "No mode")
    #     # if none
    #     else:
    #         print(sequence_list[i], ':', line)
    print('DONE\n\n')
