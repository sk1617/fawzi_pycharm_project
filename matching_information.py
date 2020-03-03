print('hi')

import glob
import pandas as pd
import numpy as np
import statistics as st
import matplotlib.pyplot as plt
import re
from DataProcessing import *
from DataImporting import sequence_list
peak_list, residue_list = main_data_processing()


def sorted_list(name):
    assert type(name) == str
    ls = list()
    for i in peak_list:
        shift = i.get_data(name)
        if type(shift) is float:
            ls.append(shift)
    ls.sort()
    return ls


def sorted_delta_list(ls):
    ls1 = [round(min([abs(x - shift) for x in ls if x != shift]), 4) for shift in ls]
    ls1.sort()
    return ls1


def find_closest(shift_list, shift_prime_list):
    ls = []
    for shift in shift_list:
        ls1 = [abs(shift - prime) for prime in shift_prime_list]
        ls.append(min(ls1))
    return ls

# print(sorted_delta_list(sorted_list('CAShift')))
# print(sorted_delta_list(sorted_list('CBShift')))

# CA Shift
x = find_closest(sorted_list('CAShift'), sorted_list('CAPrimeShift'))
plt.hist(x, bins=20)
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, max(x), step=.05))
plt.ylabel('Count')
plt.title('CAShift Deltas')
plt.show()

# CA Truncated
z = find_closest(sorted_list('CAShift'), sorted_list('CAPrimeShift'))
x = [y for y in z if y < .5]
plt.hist(x, bins=20, color='b')
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, .5, step=.05))
plt.ylabel('Count')
plt.title('TRUNCATED CAShift Deltas')
plt.show()

# CB Shift
x = find_closest(sorted_list('CBShift'), sorted_list('CBPrimeShift'))
plt.hist(x, bins=20, color='r')
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, max(x), step=.2))
plt.ylabel('Count')
plt.title('CBShift Deltas')
plt.show()

# CB Shift TRUNCATED
z = find_closest(sorted_list('CBShift'), sorted_list('CBPrimeShift'))
x = [y for y in z if y < .5]

plt.hist(x, bins=20, color='m')
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, max(x)+.05, step=.05))
plt.ylabel('Count')
plt.title('TRUNCATED CBShift Deltas')
plt.show()

# NOESY List
NOESY_shift_list = []
for peak in peak_list:
    [NOESY_shift_list.append(shift) for shift in peak.get_data('NOESYHShift')]
print(NOESY_shift_list[0:100])

# H Shift
x = find_closest(sorted_list('TROSYHShift'), NOESY_shift_list)
plt.hist(x, bins=20, color='c')
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, max(x), step=.005))
plt.ylabel('Count')
plt.title('H Shift (NOESY) Deltas')
plt.show()
print(max(x))

# H Shift TRUNCATED
x = [y for y in sorted_delta_list(sorted_list('TROSYHShift')) if y < .1]
plt.hist(x, bins=20, color='y')
plt.xlabel('Closest Shift Delta (ppm)')
plt.xticks(np.arange(0, max(x), step=.01))
plt.ylabel('Count')
plt.title('H Shift (NOESY) Deltas')
plt.show()

########################################################################################################################

two_d_dict={}
two_d_dict['HShift'] = [x.get_data('TROSYHShift') for x in peak_list]
two_d_dict['NShift'] = [x.get_data('TROSYNShift') for x in peak_list]
two_d_dict['CA'] = [x.get_data('CAShift') for x in peak_list]
two_d_dict['CB'] = [x.get_data('CBShift') for x in peak_list]
two_d_dict['CAPrime'] = [x.get_data('CAPrimeShift') for x in peak_list]
two_d_dict['CBPrime'] = [x.get_data('CBPrimeShift') for x in peak_list]
two_d_dataframe = pd.DataFrame(two_d_dict)
print(two_d_dataframe.head())

closest_ca_list, closest_cb_list = [], []
for row1 in range(len(two_d_dataframe)):
    ca_shift = two_d_dataframe.loc[row1, "CA"]
    cb_shift = two_d_dataframe.loc[row1, "CB"]
    if np.isnan(ca_shift):
        closest_ca_list.append(None)
    else:
        ls_ca = []
        for row2 in range(len(two_d_dataframe)):
            ca_prime_shift = two_d_dataframe.loc[row2, "CAPrime"]
            if not np.isnan(ca_prime_shift):
                ls_ca.append(round(abs(ca_shift - ca_prime_shift), 4))
        ls_ca.sort()
        closest_ca_list.append(ls_ca[0:5])

    if np.isnan(cb_shift):
        closest_cb_list.append(None)
    else:
        ls_cb = []
        for row2 in range(len(two_d_dataframe)):
            cb_prime_shift = two_d_dataframe.loc[row2, "CBPrime"]
            if not np.isnan(cb_prime_shift):
                ls_cb.append(round(abs(cb_shift - cb_prime_shift), 4))
        ls_cb.sort()
        closest_cb_list.append(ls_cb[0:5])

# print(closest_ca_list[0:3])
# print(closest_cb_list[0:3])
two_d_dataframe['ClosestCAShift1D'] = closest_ca_list
two_d_dataframe['ClosestCBShift1D'] = closest_cb_list

closest_2d_list = []
for row1 in range(len(two_d_dataframe)):
    ca_shift = two_d_dataframe.loc[row1, "CA"]
    cb_shift = two_d_dataframe.loc[row1, "CB"]
    if np.isnan(ca_shift):
        closest_2d_list.append(None)
    elif np.isnan(cb_shift):
        closest_2d_list.append(None)
    else:
        ls = []
        for row2 in range(len(two_d_dataframe)):
            ca_prime_shift = two_d_dataframe.loc[row2, "CAPrime"]
            cb_prime_shift = two_d_dataframe.loc[row2, "CBPrime"]
            if not np.isnan(ca_prime_shift) and not np.isnan(cb_prime_shift):
                ca_delta_squared = (abs(ca_shift - ca_prime_shift)) ** 2
                cb_delta_squared = (abs(cb_shift - cb_prime_shift)) ** 2
                ls.append((ca_delta_squared + cb_delta_squared) ** .5)
        ls.sort()
        closest_2d_list.append([round(x, 5) for x in ls[0:5]])
two_d_dataframe['Closest2DShift (ppm)'] = closest_2d_list

