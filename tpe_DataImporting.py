sequence_list = 'MSEYIRVTEDENDE'
import random

num_to_return = .1
def returns_number():
    return random.random()


energy_if_false = 200.

aa_p_d_delta = .05
aa_p_d_sn_factor = .186
aa_p_d_sub = 0
aa_p_d = 30

bb_p_d_delta = .05
bb_p_d_sn_factor = .316
bb_p_d_sub = 0
bb_p_d = 10

bmrb = 1.
a_bmrb_sn_factor = .186
b_bmrb_sn_factor = .316

noesy_perfect_match_threshold = .02
npmt_penatly = 0.0
noesy_semi_perfect_match_threshold = .05
nspmt_penalty= 400.
n_no_match_penalty = 900.
# @profile
def dist_factor_forumla(dist):
    if dist < 3:
        return 2.8
    elif dist >= 3 and dist < 7:
        return 25/(dist ** 2)
    elif dist >= 7:
        return .5
# @profile
def sn_factor_formula(sn):
    return .22 * (sn ** .5)


iterations = 1e6
temp = 1e6
exponent = .7


iterations = int(iterations)
should_append_DL = True if iterations <= 1e3 else False

# NMR experiments
if True:
    trosy = open('Test 2n4p peak list/TDP1-80-S48E-HSQC_tpe.txt')
    hnca = open('Test 2n4p peak list/TDP1-80-S48E-HNCA_tpe.txt')
    hncoca = open('Test 2n4p peak list/TDP1-80-S48E-HNCOCA_tpe.txt')
    hncocacb = open('Test 2n4p peak list/TDP1-80-S48E-HNCOCACB_tpe.txt')
    hncacb = open('Test 2n4p peak list/TDP1-80-S48E-HNCACB_tpe.txt')
    noesy = open('Test 2n4p peak list/TDP1-80-S48E-NOESY_tpe.txt')

    ftolist_file = open('Test 2n4p peak list/2n4p distances_tpe.txt')

    TROSY = list()
    for line in trosy.readlines():
        i = line.split()
        i.insert(3, 1.)
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        addition = tuple(i[:])
        TROSY.append(addition)

    HNCA = list()
    for line in hnca.readlines():
        i = line.split()
        i.insert(4, 1.)
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCA.append(addition)

    HNCOCA = list()
    for line in hncoca.readlines():
        i = line.split()
        i.insert(4, 1.)
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCOCA.append(addition)

    HNCOCACB = list()
    for line in hncocacb.readlines():
        i = line.split()
        i.insert(4, 1.)
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCOCACB.append(addition)

    HNCACB = list()
    for line in hncacb.readlines():
        i = line.split()
        i.insert(4, 1.)
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCACB.append(addition)

    # HNCO = list()
    # for line in hnco.readlines():
    #     i = line.split()
    #     i[1] = float(i[1])
    #     i[2] = float(i[2])
    #     i[3] = float(i[3])
    #     i[4] = float(i[4])
    #     i[5] = float(i[5])
    #     addition = tuple(i[:])
    #     HNCO.append(addition)

    NOESY = list()
    for line in noesy.readlines():
        i = line.split()
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        NOESY.append(addition)

    FTOList = list()
    for line in ftolist_file.readlines():
        split = line.split()
        split[0], split[1] = split[1], split[0]
        split[3], split[4] = split[4], split[3]
        split[0] = int(split[0])
        split[3] = int(split[3])
        split[6] = float(split[6])
        FTOList.append(split)


