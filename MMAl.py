from EnergyFunctions import *
from DataProcessing import *
from DataImporting import *
import random
import math
import time
import datetime as dt

Delta_List = []; energy_list = []
# @profile
def main(iterations, temp, exponent, energy_if_false, peak_list, residue_list, create_chain = False):
    start_time = time.time()  # notes start time

    print('iterations:', iterations, 'temp:', temp, 'exponent:', exponent, 'eif:', energy_if_false)  # prints initial
    if should_save_to_file:
        fh.write('iterations:{}, temp:{}, exponent:{}, eif:{}\n'.format(iterations, temp, exponent, energy_if_false))

    index_list = list(range(len(peak_list)))  # creates index list to link residue and peak lists
    assert len(index_list) == len(residue_list)
    random.shuffle(index_list, lambda: random.random())  # shuffles index list
    og_index_list = index_list  # stores original index list
    count_temp = 0  # count temp is reset
    count = 0  # counts iterations

    # iterations
    for i in range(iterations):
        count += 1
        count_temp += 1
        a = random.randint(0, len(index_list) - 1)  # select a random place from the index list
        new_index_list = index_list[:]  # makes a deep copy

        if create_chain:
            r, f = chain_creator(index_list, peak_list, a)  # gets how big the chain should be

            if True:
                a_list = [x for x in range(a + r, a + f + 1)]  # all 'a' we're switching
                while True:
                    b = random.randint(-1 * r, len(index_list) - 1 - f)  # returns 'b' inclusive
                    b_list = [x for x in range(b + r, b + f + 1)]
                    if set(a_list).isdisjoint(b_list):  # Return True if two sets have a null intersection.
                        break
                new_index_list[a + r: a + f + 1], new_index_list[b + r: b + f + 1], = \
                    new_index_list[b + r: b + f + 1], new_index_list[a + r: a + f + 1]  # swaps

        else:
            b = random.randint(0, len(index_list) - 1)
            new_index_list[a], new_index_list[b] = new_index_list[b], new_index_list[a]
            a_list = [a]; b_list = [b]

        if swap_evaluator(a_list, b_list, index_list, new_index_list, temp,
                          peak_list, residue_list, energy_if_false, Delta_List):
            index_list = new_index_list

        if temp == 0:
            print('hi')

        # TODO: FIX THIS SO IT DOESNT TAKE UP TIME
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

        print('final temperature: {}'.format(temp))
        if should_save_to_file: fh.write('final temperature: {}'.format(temp))

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
def chain_creator(index_list, peak_list, a):
    # 120 is probably a good cutoff
    f = 0  # forward
    r = 0  # reverse
    if True:
        while True:
            i = a + f  # recall, a is the random position chosen on the index list, so i is the most forward index
            if i >= len(index_list) - 1:  # ensures chain doesn't extend beyond index list
                break
            ca_ca_prime_diff = CA_CA_prime_diff(i, index_list[i], index_list, peak_list, [])  # gets energy for next pk
            cb_cb_prime_diff = CB_CB_prime_diff(i, index_list[i], index_list, peak_list, [])
            sum_energy = ca_ca_prime_diff + cb_cb_prime_diff
            if sum_energy < 600:
                p = math.exp(-0.007 * sum_energy)
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
            sum_energy = ca_ca_prime_diff + cb_cb_prime_diff
            if sum_energy < 600:
                p = math.exp(-0.007 * sum_energy)
                if random.random() < p:
                    r -= 1
                else:
                    break
            else:
                break

    assert r <= 0; assert f >= 0
    return r, f


# HOW LONG DOES THE PROGRAM RUN??? (you can also add a specific time instead of a delta)
date = dt.datetime.now() + dt.timedelta(days=6, hours=23, minutes=50)
if __name__ == '__main__':
    if dt.datetime.now() > date:
        raise Exception('end time before begin time')

    should_save_to_file = True  # SHOULD IT WRITE TO A FILE
    if should_save_to_file:
        file_name = str(dt.datetime.now().isoformat())
        fh = open('mmals/' + file_name + '.txt', 'w+')

    # starts program
    if True:
        time_start = dt.datetime.now()
        print('program started at:{}'.format(time_start))
        peak_list, residue_list = main_data_processing()
        for i in range(2):  # HOW MANY TRIALS SHOULD HAPPEN
            result = main(int(iterations), temp, exponent, energy_if_false, peak_list, residue_list, create_chain=True)

            print('Done\n\n')
            if should_save_to_file:
                fh.write('Done\n\n')

            if dt.datetime.now() > date:
                break
        time_finish = dt.datetime.now()
        print('program ended at: {} \n time delta: {}'.format(time_finish, time_finish - time_start))

    if should_save_to_file:
        fh.close()

