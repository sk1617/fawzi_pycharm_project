from MMAl import *
from EnergyFunctions import *
from tpe_DataProcessing import main_data_processing
import matplotlib.pyplot as plt
import statistics as st
import math
import pandas as pd

peak_list, residue_list = main_data_processing()
perfect_list = [x for x in range(74)]

it = sequence_list.__iter__()
w = 80
v = 0
for u, thing in enumerate(it):
    if thing == 'P':
        v += 1
        perfect_list.insert(u, w - v)


def prob(ca_ca_prime_diff, cb_cb_prime_diff):
    if ca_ca_prime_diff and cb_cb_prime_diff and \
            ca_ca_prime_diff < 900 and cb_cb_prime_diff < 900:
        p = ((1 + math.exp(.01 * (ca_ca_prime_diff - 500))) ** -1) * \
            ((1 + math.exp(.04 * (cb_cb_prime_diff - 400))) ** -1)
        return p
    else:
        return False

prob_list = []
for i, pai in enumerate(perfect_list):
    ca_ca_prime_diff = CA_CA_prime_diff(i, pai, perfect_list, peak_list, [])
    cb_cb_prime_diff = CB_CB_prime_diff(i, pai, perfect_list, peak_list, [])
    p = prob(ca_ca_prime_diff, cb_cb_prime_diff)
    prob_list.append((pai, residue_list[i].get_data('residue_name'), p,
                      ca_ca_prime_diff, cb_cb_prime_diff))

df = pd.DataFrame(prob_list, columns = ['pai', 'name', 'p', 'ca', 'cb'])
# df


writer = pd.ExcelWriter('chain creator for p 1.xlsx')
df.to_excel(writer, 'perfect list 1')
writer.close()

