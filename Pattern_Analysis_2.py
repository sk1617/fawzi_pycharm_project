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

    # index list specific
    data = pd.DataFrame(columns=['residue number', 'amino acid', 'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime',
                             'CACAPrime', 'CBCBPrime', 'CACBExp', 'CACBExpPrime', 'NOESY', 'NearbyHShiftList', 'ClosestMatchDelta'])
    if True:
        # index_list = initial_data.loc[158, 'index_list']
        index_list = [378, 415, 139, 313, 8, 423, 273, 380, 266, 26, 60, 356, 73, 94, 322, 329, 245, 330, 53, 292, 38, 77, 290, 466, 134, 426, 155, 176, 64, 197, 83, 232, 431, 30, 173, 105, 11, 317, 75, 136, 120, 425, 316, 304, 15, 324, 72, 18, 210, 55, 448, 182, 472, 113, 343, 405, 167, 115, 369, 85, 470, 187, 473, 152, 172, 132, 474, 52, 463, 203, 235, 411, 170, 62, 364, 388, 125, 141, 437, 9, 81, 208, 278, 410, 309, 189, 31, 166, 300, 348, 464, 218, 181, 310, 376, 303, 216, 101, 2, 194, 192, 236, 14, 264, 193, 392, 243, 336, 428, 390, 63, 446, 98, 142, 337, 128, 237, 467, 289, 421, 250, 106, 200, 374, 409, 151, 379, 460, 3, 301, 112, 28, 44, 359, 287, 270, 100, 217, 417, 308, 434, 12, 13, 358, 256, 347, 457, 46, 293, 183, 267, 164, 274, 366, 48, 401, 178, 280, 45, 430, 286, 323, 325, 20, 169, 10, 114, 196, 260, 252, 351, 449, 99, 439, 135, 298, 71, 328, 168, 413, 191, 299, 258, 349, 296, 370, 90, 327, 408, 234, 195, 344, 263, 455, 357, 137, 97, 244, 24, 346, 221, 209, 174, 282, 246, 288, 468, 36, 402, 454, 57, 342, 41, 407, 230, 386, 432, 156, 111, 88, 318, 338, 131, 87, 350, 33, 214, 272, 69, 238, 110, 381, 305, 47, 143, 360, 334, 49, 435, 403, 283, 452, 16, 418, 404, 70, 340, 147, 22, 394, 51, 312, 186, 219, 456, 163, 148, 161, 419, 43, 84, 175, 7, 92, 80, 451, 353, 79, 306, 382, 154, 65, 377, 383, 124, 146, 220, 91, 459, 190, 61, 0, 160, 319, 367, 34, 389, 375, 126, 171, 66, 93, 180, 162, 269, 332, 284, 138, 198, 395, 339, 253, 223, 56, 429, 205, 233, 25, 458, 372, 424, 121, 321, 257, 400, 226, 150, 261, 239, 393, 387, 249, 345, 368, 450, 123, 133, 157, 59, 206, 117, 129, 29, 362, 465, 222, 251, 442, 188, 453, 231, 103, 384, 179, 241, 385, 433, 225, 441, 259, 229, 291, 444, 184, 438, 331, 27, 297, 104, 185, 130, 371, 1, 165, 277, 265, 23, 116, 74, 37, 42, 295, 127, 365, 58, 406, 76, 302, 68, 177, 352, 275, 422, 5, 122, 320, 268, 307, 89, 461, 361, 6, 436, 227, 355, 262, 427, 363, 255, 445, 228, 414, 109, 471, 285, 420, 341, 78, 145, 213, 314, 281, 159, 17, 204, 271, 224, 82, 144, 19, 50, 416, 202, 447, 40, 207, 119, 54, 242, 86, 201, 440, 373, 294, 211, 248, 240, 254, 469, 102, 443, 397, 396, 391, 212, 276, 311, 279, 35, 333, 32, 118, 158, 107, 96, 21, 95, 412, 140, 462, 315, 39, 149, 399, 4, 199, 108, 398, 335, 354, 67, 153, 326, 247, 215]
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
            nearby_h_shift, noesy_delta = [], []

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
                                nearby_h_shift = info_list[k + l][2]
                                noesy_delta = info_list[k + l][3]
                                noesy_delta = [round(x, 4) for x in noesy_delta]

            temp_df = pd.DataFrame([[residue_number, amino_acid, n_shift, h_shift, ca_shift, ca_prime_shift, cb_shift, cb_prime_shift,
                                     ca_ca_prime, cb_cb_prime, ca_cb_exp, ca_cb_exp_prime, noesy_energy, nearby_h_shift, noesy_delta]],
                                   columns=['residue number', 'amino acid', 'N', 'H', 'CA', 'CAPrime', 'CB', 'CBPrime',
                                            'CACAPrime', 'CBCBPrime', 'CACBExp', 'CACBExpPrime', 'NOESY', 'NearbyHShiftList', 'ClosestMatchDelta'])
            data = data.append(temp_df, ignore_index=True)
        print(data.head())
        data.to_csv('/Volumes/Transcend/fawzi_pycharm_project/data_TEST_{}.csv'.format('2021-01-03'))
