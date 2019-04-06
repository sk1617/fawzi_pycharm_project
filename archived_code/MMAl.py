from archived_code.EnergyFunctions import *
from DataProcessing import *
from DataImporting import *
import random
import math


def main(iterations, temp, exponent):
    # makes index list
    index_list = []
    peak_and_residue_list = main_data_processing()
    peak_list, residue_list = peak_and_residue_list[0], peak_and_residue_list[1]
    randomized_peak_list = peak_list
    random.shuffle(randomized_peak_list, returns_number)
    for i, residue in enumerate(residue_list):
        index_list.append([randomized_peak_list[i], residue])

    og_index_list = index_list
    count_temp = 0
    count = 0

    for i in range(iterations):
        if temp == 0.0:
            break
        count += 1
        count_temp += 1
        a = random.randint(0, len(index_list) - 1)
        b = random.randint(0, len(index_list) - 1)
        new_index_list = index_list[:]
        new_index_list[a][0], new_index_list[b][0] = new_index_list[b][0], new_index_list[a][0]

        if swap_evaluator(index_list, new_index_list, temp):
            index_list = new_index_list

        if count_temp == int(iterations*((math.log(1/temp, exponent))**-1)):
            count_temp = 0
            temp = temp * exponent

    return(og_index_list, index_list)



def swap_evaluator(index_list, new_index_list, temp):
    new_energy, old_energy = float(), float()

    delta = main_energy_function(a, b, old_index_list, new_index_list, peak_list, residue_list, energy_if_false, Delta_List)

    if delta <= 0:
        return True
    elif delta > 0:
        p = math.exp(-delta/temp)
        return True if random.random() <= p else False


if __name__ == '__main__':
    pair = main(iterations, temp, exponent)