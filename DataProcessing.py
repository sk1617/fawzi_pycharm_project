from itertools import chain
import random
from EnergyFunctions import *
from DataImporting import *
import numpy as np
import statistics as st

def main_data_processing():
    peak_list = []
    create_peaks(peak_list)

    residue_list = []
    create_amino_acids(residue_list)

    return peak_list, residue_list


class Peaks:  # class and required defs for each peak
    def __init__(self, **kwargs):
        self.peak_properties = kwargs

    def add_data(self, label, value):
        self.peak_properties[label] = value

    def get_data(self, label):
        return self.peak_properties[label]

    def add_data_to_list(self, list_name, value):
        self.peak_properties[list_name].append(value)


def create_peaks(peak_list):  # instantiate every peak as a class and append it to list peak_list
    for index in range(len(sequence_list)):
        try:
            peak_list.append(Peaks(peakNumber=index,
                                   TROSYHShift=outside_peaks_df.loc[index, 'H_'],
                                   TROSYNShift=outside_peaks_df.loc[index, 'N_'],
                                   CAShift=outside_peaks_df.loc[index, 'CA'],
                                   CAPrimeShift=outside_peaks_df.loc[index, 'CA-1'],
                                   CBShift=outside_peaks_df.loc[index, 'CB'],
                                   CBPrimeShift=outside_peaks_df.loc[index, 'CB-1'],
                                   COPrimeShift=outside_peaks_df.loc[index, 'CO-1'],
                                   NOESYHShift=[], NOESYNShift=[], NearbyHShift=[]))
        except KeyError:
            peak_list.append(Peaks(peakNumber=index,
                                   TROSYHShift=None, TROSYNShift=None, CAShift=None, CAPrimeShift=None,
                                   CBShift=None, CBPrimeShift=None, COPrimeShift=None,
                                   NOESYHShift=[], NOESYNShift=[], NearbyHShift=[]))

    fix_none(peak_list)
    old_data(peak_list)


def old_data(peak_list):
    old_peak_list = []
    for i, line in enumerate(TROSY):
        old_peak_list.append(Peaks(peakNumber=i,
                                   TROSYHShift=TROSY[i][2], TROSYNShift=TROSY[i][1],
                                   HNCAHShift=None, HNCANShift=None, CAShift=None,
                                   HNCOCAHShift=None, HNCOCANShift=None, CAPrimeShift=None,
                                   HNCACBHShift=None, HNCACBNShift=None, CBShift=None,
                                   HNCOHShift=None, HNCONShift=None, COPrimeShift=None,
                                   NOESYHShift=[], NOESYNShift=[], NearbyHShift=[], NOESYSignalNoise=[]))

    for line in HNCA:
        if 'PRIME' in line[0]:
            continue
        else:
            # memory for HCNA line:
            peak_number = 0
            H, N = line[3], line[2]
            for i, peak in enumerate(old_peak_list):
                if peak.get_data('TROSYHShift') is None:
                    continue
                else:
                    prev_H = old_peak_list[peak_number].get_data('TROSYHShift')
                    prev_N = old_peak_list[peak_number].get_data('TROSYNShift')
                    new_H = peak.get_data('TROSYHShift')
                    new_N = peak.get_data('TROSYNShift')
                    if distance_formula(H, new_H, N, new_N) < distance_formula(H, prev_H, N, prev_N) and \
                            distance_formula(H, new_H, N, new_N) < .075:
                        peak_number = i

            old_peak_list[peak_number].add_data('HNCAHShift', H)
            old_peak_list[peak_number].add_data('HNCANShift', N)
            old_peak_list[peak_number].add_data('CAShift', line[1])

    for line in HNCOCA:
        peak_number = 0  # memory for HCNA line:
        H, N = line[3], line[2]
        for i, peak in enumerate(old_peak_list):
            if peak.get_data('TROSYHShift') is None:
                continue
            else:
                prev_H = old_peak_list[peak_number].get_data('TROSYHShift')
                prev_N = old_peak_list[peak_number].get_data('TROSYNShift')
                new_H = peak.get_data('TROSYHShift')
                new_N = peak.get_data('TROSYNShift')
                if distance_formula(H, new_H, N, new_N) < distance_formula(H, prev_H, N, prev_N) and \
                        distance_formula(H, new_H, N, new_N) < .075:
                    peak_number = i

        old_peak_list[peak_number].add_data('HNCOCAHShift', H)
        old_peak_list[peak_number].add_data('HNCOCANShift', N)
        old_peak_list[peak_number].add_data('CAPrimeShift', line[1])

    for line in HNCO:
        # memory for hnco line:
        peak_number = 0
        H, N = line[3], line[2]
        for i, peak in enumerate(old_peak_list):
            if peak.get_data('TROSYHShift') is None:
                continue
            else:
                prev_H = old_peak_list[peak_number].get_data('TROSYHShift')
                prev_N = old_peak_list[peak_number].get_data('TROSYNShift')
                new_H = peak.get_data('TROSYHShift')
                new_N = peak.get_data('TROSYNShift')
                if distance_formula(H, new_H, N, new_N) < distance_formula(H, prev_H, N, prev_N) and \
                        distance_formula(H, new_H, N, new_N) < .075:
                    peak_number = i

        old_peak_list[peak_number].add_data('HNCOHShift', H)
        old_peak_list[peak_number].add_data('HNCONShift', N)
        old_peak_list[peak_number].add_data('COPrimeShift', line[1])

    for line in NOESY:
        peak_number = 0
        H = line[3]
        N = line[2]
        for i, peak in enumerate(old_peak_list):  # if new is less than old
            if peak.get_data('TROSYHShift') is None:
                continue
            else:
                prev_H = old_peak_list[peak_number].get_data('TROSYHShift')
                prev_N = old_peak_list[peak_number].get_data('TROSYNShift')
                new_H = peak.get_data('TROSYHShift')
                new_N = peak.get_data('TROSYNShift')
                if distance_formula(H, new_H, N, new_N) < distance_formula(H, prev_H, N, prev_N):
                    peak_number = i

        old_peak_list[peak_number].add_data_to_list('NOESYHShift', H)
        old_peak_list[peak_number].add_data_to_list('NOESYNShift', N)
        old_peak_list[peak_number].add_data_to_list('NearbyHShift', line[1])

    for old_peak in old_peak_list:
        if len(old_peak.get_data('NOESYNShift')) > 0:
            old_h = st.mean(old_peak.get_data('NOESYHShift'))
            old_n = st.mean(old_peak.get_data('NOESYNShift'))
        else:
            old_n = old_peak.get_data('TROSYNShift')
            old_h = old_peak.get_data('TROSYHShift')

        old_co_prime = old_peak.get_data('COPrimeShift')
        memory = 10000
        distance = 10000
        if old_co_prime is None:
            continue
        for i, new_peak in enumerate(peak_list):
            new_n = new_peak.get_data('TROSYNShift')
            new_h = new_peak.get_data('TROSYHShift')
            new_co_prime = new_peak.get_data('COPrimeShift')
            if new_h is None or np.isnan(new_h):
                continue
            # elif new_co_prime is None or np.isnan(new_co_prime):
            #     continue
            else:
                newdist = distance_formula(old_h, new_h, old_n, new_n, old_co_prime, new_co_prime)
                if newdist < distance:
                    distance = newdist
                    memory = i

        new_peak = peak_list[memory]
        if abs(old_h - new_peak.get_data('TROSYHShift')) < .1:
            if abs(old_n - new_peak.get_data('TROSYNShift')) < 1:
                # NOESY
                new_peak.add_data('NOESYHShift', old_peak.get_data('NOESYHShift'))
                new_peak.add_data('NOESYNShift', old_peak.get_data('NOESYNShift'))
                new_peak.add_data('NearbyHShift', old_peak.get_data('NearbyHShift'))
                # CA Shift
                if new_peak.get_data('CAShift') is None or np.isnan(new_peak.get_data('CAShift')):
                    new_peak.add_data('CAShift', old_peak.get_data('CAShift'))
                # CA Prime Shift
                if new_peak.get_data('CAPrimeShift') is None or np.isnan(new_peak.get_data('CAPrimeShift')):
                    new_peak.add_data('CAPrimeShift', old_peak.get_data('CAPrimeShift'))


def fix_none(peak_list):
    for peak in peak_list:
        for prop in peak.peak_properties:
            if 'NOESY' in prop or 'Nearby' in prop:
                continue
            elif peak.get_data(prop) is None:
                continue
            elif np.isnan(peak.get_data(prop)):
                peak.add_data(prop, None)

'''
----------------------------------------//----------------------------------------
Now to processing sequence data
'''


class AminoAcid:  # class and required defs for each amino acid
    def __init__(self, **kwargs):
        self.residue_properties = kwargs

    def add_data(self, label, value):
        self.residue_properties[label] = value

    def get_data(self, label):
        """

        :rtype:
        """
        return self.residue_properties[label]

    def add_data_to_list(self, list_name, value):
        self.residue_properties[list_name].append(value)


def create_amino_acids(residue_list):
    for i, residue in enumerate(sequence_list):
        residue_list.append(AminoAcid(residue_name=residue, index=i + 1,
                                      secondary_structure=None,
                                      CAExpected=None, CAExpectedSD=None, CBExpected=None, CBExpectedSD=None,
                                      closebyAminoAcids=[], closebyAminoAcidsDist=[]))

    secondary_structure(residue_list)
    expected_CACB(residue_list)
    HN_distances(residue_list)


def secondary_structure(residue_list):  # defines secondary structure for all positions
    for residue in residue_list:
        # these all need to be shifted by one bc we have an additional starting serine
        random_range = [x + 1 for x in chain(range(0,6), range(15,18), range (22,22), range (26,28), range (44,64),
                                         range(70,74), range(78,80), range(84,100), range(126, 160), range(166,170),
                                         range(177,195), range(201,211), range(218,239), range(246,254), range(259,263),
                                         range(267,270),range(280,285), range(292,300), range(314,330),
                                         range (347, 347), range(354,358), range(361,366), range(391,408),
                                         range(426,461), range(465,477), range(22,23), range(347,348))]

        if (residue.get_data('index')+1) in random_range:
            residue.add_data('secondary_structure', 'R')

        # Alpha
        alpha_range = [y + 1 for y in chain(range(6,15), range(23,26), range(28,44), range(100,126),range(160,166),
                                        range(270,274), range(300,314), range(330,347), range(348,354),
                                        range(358,361), range(366,391), range(408,426), range(461,465))]

        if (residue.get_data('index')+1) in alpha_range:
            residue.add_data('secondary_structure', 'A')

        # Beta
        beta_range = [z + 1 for z in chain(range(18,22), range(64,70), range(74,78), range(80,84),
                                       range(170,177), range(195,201), range(211,218), range(239,247), range(254,259),
                                       range(263,267), range(274,280), range(285,292))]

        if (residue.get_data('index')+1) in beta_range:
            residue.add_data('secondary_structure', 'B')


def expected_CACB(residue_list):  # defines parameters of BMRB data in class BMRB
    from BMRB import A as BMRB
    """
    (ResName, CA Beta, SD, CA Random, SD, CA Alpha, SD, CB Beta, SD, CB Random, SD, CB Alpha, SD)
       0        1      2       3      4      5      6       7    8       9      10      11    12
    Values range from 0-12
    """
    for residue in residue_list:
        # assign an index based on the secondary structure. This index is used to get the right CACB shifts
        if residue.get_data('secondary_structure') == 'R':
            index = 1
        elif residue.get_data('secondary_structure') == 'A':
            index = 2
        elif residue.get_data('secondary_structure') == 'B':
            index = 0
        elif residue.get_data('secondary_structure') is None:
            continue

        for letter in range(0, 20):  # sees if the amino acid of AminoAcid matches the amino acid in the table
            if BMRB[letter][0] == residue.get_data('residue_name'):
                residue.add_data('CAExpected', BMRB[letter][1 + index * 2])
                residue.add_data('CAExpectedSD', BMRB[letter][2 + index * 2])

                residue.add_data('CBExpected', BMRB[letter][7 + index * 2])
                residue.add_data('CBExpectedSD', BMRB[letter][8 + index * 2])


def HN_distances(residue_list):
    """
    NOTE:
    FTOList is a list of lists: the first 2 elements are:
    [['1', 'THR', 'HN', '3', 'LYS+', 'HN', '4.99'], ['1', 'THR', 'HN', '4', 'ASP', 'HN', '4.61'],  ...]
    """

    for line in FTOList:
        # these all need to be shifted by one bc we have an additional starting serine
        res_1_num = int(line[0]) + 1
        res_2_num = int(line[3]) + 1
        dist = float(line[6])
        if dist < 8: #adding all residues within 8 angstroms
            # adding the line to each of the Amino Acids
            res_1 = residue_list[res_1_num - 1]  # subtract one bc it starts counting the AAL at 0
            res_2 = residue_list[res_2_num - 1]

            if res_2_num not in res_1.get_data('closebyAminoAcids'):
                res_1.add_data_to_list('closebyAminoAcids', res_2_num)
                res_1.add_data_to_list('closebyAminoAcidsDist', dist)

            if res_1_num not in res_2.get_data('closebyAminoAcids'):
                res_2.add_data_to_list('closebyAminoAcids', res_1_num)
                res_2.add_data_to_list('closebyAminoAcidsDist', dist)


def distance_formula(H1, H2, N1, N2, CO1=None, CO2=None):
    if CO1 is not None and CO2 is not None:
        z = ((((H1-H2)/8.0561) ** 2 + ((N1-N2)/120.2105) ** 2 + ((CO1-CO2)/176) ** 2) ** .5) * 2/3
    else:
        z = (((H1-H2)/8.0561) ** 2 + ((N1-N2)/120.2105) ** 2) ** .5
    return z


if __name__ == '__main__':
    pair = main_data_processing()
    peak_list, residue_list = pair[0], pair[1]
