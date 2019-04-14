from tpe_EnergyFunctions import *
from DataProcessing import *
from full_DataImporting import *
import random
import math
import time
import datetime as dt

Delta_List = []; energy_list = []
# @profile
def main(iterations, temp, exponent, energy_if_false, peak_list, residue_list, create_chain = False):
    start_time = time.time()
    if dt.datetime.now() < date:
        print('iterations:', iterations, 'temp:', temp, 'exponent:', exponent, 'eif:', energy_if_false)
        if should_save_to_file:
            fh.write('iterations:{}, temp:{}, exponent:{}, eif:{}\n'.format(iterations, temp, exponent, energy_if_false))

    index_list = list(range(len(peak_list)))
    assert len(index_list) == len(residue_list)
    random.shuffle(index_list, returns_number)

    og_index_list = index_list
    count_temp = 0
    count = 0

    # iterations
    for i in range(iterations):
        count += 1
        count_temp += 1
        a = random.randint(0, len(index_list) - 1)
        new_index_list = index_list[:]
        should_swap_eval = True

        if create_chain:
            r, f = chain_creator(index_list, peak_list, a)

            if True:
                a_list = [x for x in range(a + r, a + f + 1)]
                while True:
                    b = random.randint(-1 * r, len(index_list) - 1 - f)
                    b_list = [x for x in range(b + r, b + f + 1)]
                    if set(a_list).isdisjoint(b_list):
                        break
                new_index_list[a + r: a + f + 1], new_index_list[b + r: b + f + 1], = \
                    new_index_list[b + r: b + f + 1], new_index_list[a + r: a + f + 1]
            # set
            if True:
                if len(set(new_index_list)) != len(new_index_list):
                    print('ERROR!!')
                    print('a list: {} \n b list: {} \n r: {} \n f: {} \n index list {} \n new index list {}'.format(
                        a_list, b_list, r, f, index_list, new_index_list))
                    should_swap_eval = False

            # list
            if True:
                a_temporary = a_list[:-1]
                b_temporary = b_list[:-1]
                if not set(a_temporary).isdisjoint(b_temporary):
                    print('ERROR!!')
                    print('a list: {} \n b list: {} \n r: {} \n f: {} \n index list {} \n new index list {}'.format(
                        a_list, b_list, r, f, index_list, new_index_list))
                    should_swap_eval = False

        else:
            b = random.randint(0, len(index_list) - 1)
            new_index_list[a], new_index_list[b] = new_index_list[b], new_index_list[a]
            a_list = [a]; b_list = [b]

        if should_swap_eval:
            if swap_evaluator(a_list, b_list, index_list, new_index_list, temp, peak_list, residue_list, energy_if_false, Delta_List):
                index_list = new_index_list

            if temp == 0:
                print('hi')
            #FIX THIS SO IT DOESNT TAKE UP TIME
            elif count_temp == int(iterations*((math.log(1/temp, exponent))**-1)):
                count_temp = 0
                temp = temp * exponent

            if dt.datetime.now() > date: break

    # Data for runtime
    if True:
        end_time = time.time()
        print('time taken (sec): {}'.format(end_time - start_time))
        if should_save_to_file: fh.write('time taken (sec): {}\n'.format(end_time - start_time))

        print('number of iterations: {}'.format(count))
        if should_save_to_file: fh.write('number of iterations: {}\n'.format(count))

        print('og index list:', [x for x in og_index_list])
        if should_save_to_file: fh.write('og index list: {}\n'.format([x for x in og_index_list]))

        og_energy = eval_energy(og_index_list, peak_list, residue_list, energy_if_false)[0]
        print('og index list energy:', og_energy)
        if should_save_to_file: fh.write('og index list energy: {}\n'.format(og_energy))

        print('index list:', [x for x in index_list])
        if should_save_to_file: fh.write('index list: {}\n'.format([x for x in index_list]))

        energy = eval_energy(index_list, peak_list, residue_list, energy_if_false)[0]
        energy_list.append(energy)
        print('index list energy:', energy)
        if should_save_to_file: fh.write('index list energy: {}'.format(energy))

    # RETURNS HERE
    return og_index_list, og_energy, index_list, energy

# @profile
def swap_evaluator(a_list, b_list, index_list, new_index_list, temp, peak_list, residue_list, energy_if_false, Delta_List):
    global delta
    delta = main_energy_function(a_list, b_list, index_list, new_index_list, peak_list, residue_list, energy_if_false, Delta_List)
    if should_append_DL:
        Delta_List.append([delta, a_list, b_list, [index_list[x] for x in a_list], [index_list[x] for x in b_list], index_list, new_index_list])

    if delta <= 0:
        return True
    elif temp == 0:
        return False
    elif delta > 0:
        p = math.exp(-delta/temp)
        if random.random() <= p:
            return True
        else:
            return False

# @profile
def  chain_creator(index_list, peak_list, a):
    # 120 is probably a good cutoff
    # forwards
    f = 0
    r = 0
    if True:
        while True:
            i = a + f
            if i > len(index_list):
                break
            ca_ca_prime_diff = CA_CA_prime_diff(i, index_list[i], index_list, peak_list, [])
            cb_cb_prime_diff = CB_CB_prime_diff(i, index_list[i], index_list, peak_list, [])
            if ca_ca_prime_diff and cb_cb_prime_diff and ca_ca_prime_diff < 900 and cb_cb_prime_diff < 900:
                p = ((1 + math.exp(.01 * (ca_ca_prime_diff - 500))) ** -1) * \
                    ((1 + math.exp(.04 * (cb_cb_prime_diff - 400))) ** -1)
                if random.random() < p:
                    f += 1
                else:
                    break
            else:
                break

        while True:
            i = a + r - 1
            if i < 0:
                break
            ca_ca_prime_diff = CA_CA_prime_diff(i, index_list[i], index_list, peak_list, [])
            cb_cb_prime_diff = CB_CB_prime_diff(i, index_list[i], index_list, peak_list, [])
            if ca_ca_prime_diff and cb_cb_prime_diff and ca_ca_prime_diff < 900 and cb_cb_prime_diff < 900:
                p = ((1 + math.exp(.01 * (ca_ca_prime_diff - 500))) ** -1) * \
                    ((1 + math.exp(.04 * (cb_cb_prime_diff - 400))) ** -1)
                if random.random() < p:
                    r -= 1
                else:
                    break
            else:
                break

    assert r <= 0; assert f >= 0
    return r, f


date = dt.datetime.now() + dt.timedelta(hours=16, minutes=58)
if __name__ == '__main__':
    if dt.datetime.now() > date:
        raise Exception('end time before begin time')

    should_save_to_file = True
    if should_save_to_file:
        file_name = str(dt.datetime.now().isoformat())
        file_name = file_name[5:-5]
        fh = open('mmals/' + file_name + '.txt', 'w+')

    # starts program
    if True:
        peak_list, residue_list = main_data_processing()
        for i in range(50):
            result = main(int(iterations), temp, exponent, energy_if_false, peak_list, residue_list, create_chain=True)

            print('Done\n\n')
            if should_save_to_file:
                fh.write('Done\n\n')

            if dt.datetime.now() > date:
                break

    if should_save_to_file:
        fh.close()


