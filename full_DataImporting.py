sequence_list = 'MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVEQCMRGVRLVEGILHAPDAGWGNLVYVVNYPKD'
import random

def returns_number():
    return random.random()

energy_if_false = 300.

aa_p_d_delta = .15
aa_p_d_sn_factor = .186
aa_p_d_sub = 0
aa_p_d = 30

bb_p_d_delta = .15
bb_p_d_sn_factor = .316
bb_p_d_sub = 0
bb_p_d = 20

bmrb = 1.5
a_bmrb_sn_factor = .186
b_bmrb_sn_factor = .316

noesy_perfect_match_threshold = .02
npmt_penatly = 0.1
noesy_semi_perfect_match_threshold = .05
nspmt_penalty = 600.
n_no_match_penalty = 1200.

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


iterations = 1e3
temp = 1e6
exponent = .7


iterations = int(iterations)
should_append_DL = True if iterations <= 1e3 else False

# NMR experiments
if True:
    tr = open('Full 2np4 peak list/TDP1-80-S48E-HSQC copy.txt', 'r')
    f1 = open('Full 2np4 peak list/TDP1-80-S48E-HNCACB copy.txt', 'r')
    f2 = open('Full 2np4 peak list/TDP1-80-S48E-CBCACONH copy.txt', 'r')

    TROSY = []
    HNCA = []
    HNCACB = []
    HNCOCA = []
    HNCOCACB = []

    for line in f1.readlines():
        split = line.split()
        if len(split) == 0:
            continue
        elif split[0][-6:] == 'CB-N-H':
            split.insert(4, 1.)
            for i in range(1, 6):
                split[i] = float(split[i])
            HNCACB.append(split)
        elif split[0][-6:] == 'CA-N-H':
            split.insert(4, 1.)
            for i in range(1, 6):
                split[i] = float(split[i])
            HNCA.append(split)

    for line in f2.readlines():
        split = line.split()
        if len(split) == 0:
            continue
        elif 'CA-' in split[0]:
            split.insert(4, 1.)
            for i in range(1, 6):
                split[i] = float(split[i])
            HNCOCA.append(split)

        elif 'CB-' in split[0]:
            split.insert(4, 1.)
            for i in range(1, 6):
                split[i] = float(split[i])
            HNCOCACB.append(split)

    for line in tr.readlines():
        split = line.split()
        if len(split) == 0:
            continue
        elif split[0] == 'Assignment':
            continue
        else:
            split.insert(3, 1.)
            for i in range(1, 5):
                split[i] = float(split[i])
            TROSY.append(split)

    FTOList = list()
    ftolist_file = open('Full 2np4 peak list/2n4p distances copy.txt')
    for line in ftolist_file.readlines():
        split = line.split()
        if float(split[6]) < 6:
            split[0], split[1] = split[1], split[0]
            split[3], split[4] = split[4], split[3]
            split[0] = int(split[0])
            split[3] = int(split[3])
            split[6] = float(split[6])
            FTOList.append(split)

NOESY = []
for closeby_line in FTOList:
    pk_num = closeby_line[0]
    close_pk_num = closeby_line[3]
    for line1 in TROSY:
        x = line1[0].index('N')
        if line1[0][1:x] == str(pk_num):
            pk_H = line1[2]
            pk_N = line1[1]
            break

    for line2 in TROSY:
        x = line2[0].index('N')
        if line2[0][1:x] == str(close_pk_num):
            clpk_H = line2[2]
            break

    app = ('{}N-NH'.format(pk_num), clpk_H, pk_N, pk_H, 0.0, 20.0)
    NOESY.append(app)