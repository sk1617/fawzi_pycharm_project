# Beginning part of the sequence list: MGSSHHHHHHSSGLVPRGSHM'
sequence_list = 'TPKDDEFYQQWQLKYPKLILREASSVSEELHKEVQEAFLTLHKHGCLFRDLVRIQGKDLLTPVSRILIGNPGCTYKYLNTRLFTVPWPVKGSNIKHTEAE' \
                'IAAACETFLKLNDYLQIETIQALEELAAKEKANEDAVPLCMSADFPRVGMGSSYNGQDEVDIKSRAAYNVTLLNFMDPQKMPYLKEEPYFGMGKMAVSWH' \
                'HDENLVDRSAVAVYSYSCEGPEEESEDDSHLEGRDPDIWHVGFKISWDIETPGLAIPLHQGDCYFMLDDLNATHQHCVLAGSQPRFSSTHRVAECSTGTL' \
                'DYILQRCQLALQNVCDDVDNDDVSLKSFEPAVLKQGEEIHNEVEFEWLRQFWFQGNRYRKCTDWWCQPMAQLEALWKKMEGVTNAVLHEVKREGLPVEQR' \
                'NEILTAILASLTARQNLRREWHARCQSRIARTLPADQKPECRPYWEKDDASMPLPFDLTDIVSELRGQLLEAKP'
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


iterations = 5e9
temp = 1e6
exponent = .9

iterations = int(iterations)
should_append_DL = True if iterations <= 1e3 else False


# NMR experiments
if True:
    trosy = open('FTO_peaklists/NH_TROSY copy.list')
    hnca = open('FTO_peaklists/hnca copy.list')
    hncoca = open('FTO_peaklists/hncoca copy.list')
    hncocacb = open('FTO_peaklists/hncocacb copy.list')
    hncacb = open('FTO_peaklists/hncacb copy.list')
    hnco = open('FTO_peaklists/hnco copy.list')
    noesy = open('FTO_peaklists/HHN_NOESY_3D copy.list')

    ftolist_file = open('FTO_peaklists/HNdistances_FTO.txt')

    TROSY = list()
    for line in trosy.readlines():
        i = line.split()
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        addition = tuple(i[:])
        TROSY.append(addition)

    HNCA = list()
    for line in hnca.readlines():
        i = line.split()
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
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCACB.append(addition)

    HNCO = list()
    for line in hnco.readlines():
        i = line.split()
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[3] = float(i[3])
        i[4] = float(i[4])
        i[5] = float(i[5])
        addition = tuple(i[:])
        HNCO.append(addition)

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
        FTOList.append(split)
