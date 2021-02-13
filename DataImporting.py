import random
import pandas as pd
# Beginning part of the sequence list: MGSSHHHHHHSSGLVPRGSHM'
sequence_list = 'STPKDDEFYQQWQLKYPKLILREASSVSEELHKEVQEAFLTLHKHGCLFRDLVRIQGKDLLTPVSRILIGNPGCTYKYLNTRLFTVPWPVKGSNIKHTEAE' \
                'IAAACETFLKLNDYLQIETIQALEELAAKEKANEDAVPLCMSADFPRVGMGSSYNGQDEVDIKSRAAYNVTLLNFMDPQKMPYLKEEPYFGMGKMAVSWH' \
                'HDENLVDRSAVAVYSYSCEGPEEESEDDSHLEGRDPDIWHVGFKISWDIETPGLAIPLHQGDCYFMLDDLNATHQHCVLAGSQPRFSSTHRVAECSTGTL' \
                'DYILQRCQLALQNVCDDVDNDDVSLKSFEPAVLKQGEEIHNEVEFEWLRQFWFQGNRYRKCTDWWCQPMAQLEALWKKMEGVTNAVLHEVKREGLPVEQR' \
                'NEILTAILASLTARQNLRREWHARCQSRIARTLPADQKPECRPYWEKDDASMPLPFDLTDIVSELRGQLLEAKP'


def returns_number():
    return random.random()

energy_if_false = 300.

# These will all need to be changed now that we have real CB data
aa_p_d = 2000
bb_p_d = 2000
bmrb_ca = 3000
bmrb_cb = 3000
noesy_weight = 3000

iterations = 1e8

temp = 1e6
exponent = .9

iterations = int(iterations)
should_append_DL = True if iterations <= 1e3 else False

outside_peaks_df = pd.read_table('fto_assignment_files/Full_length/CS29_FTO',
                                 sep='\t', index_col=False)
outside_peaks_df['N_'] = [float(x)+.23 for x in outside_peaks_df['N']]
outside_peaks_df['H_'] = [float(x)+.08 for x in outside_peaks_df['H']]
outside_peaks_df['CO-1'] = [float(x) if x != '-' else None for x in outside_peaks_df['CO-1']]
outside_peaks_df['CA'] = [float(x) if x != '-' else None for x in outside_peaks_df['CA']]
outside_peaks_df['CA-1'] = [float(x) if x != '-' else None for x in outside_peaks_df['CA-1']]
outside_peaks_df['CB'] = [float(x) if x != '-' else None for x in outside_peaks_df['CB']]
outside_peaks_df['CB-1'] = [float(x) if x != '-' else None for x in outside_peaks_df['CB-1']]

# OLD NMR experiments
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