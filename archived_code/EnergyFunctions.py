from DataImporting import *
import math


def main_energy_function(i, assignment, index_list):
    energy = 0.0

    CA_CA_prime_diff_sub_energy = CA_CA_prime_diff(i, assignment, index_list)
    if CA_CA_prime_diff_sub_energy:
        energy += CA_CA_prime_diff_sub_energy
    else:
        energy += energy_if_false

    CB_CB_prime_diff_sub_energy = CB_CB_prime_diff(i, assignment, index_list)
    if CB_CB_prime_diff_sub_energy:
        energy += CB_CB_prime_diff_sub_energy
    else:
        energy += energy_if_false

    CA_obs_BMRB_diff_sub_energy = CA_obs_BMRB_diff(i, assignment)
    if CA_obs_BMRB_diff_sub_energy:
        energy += CA_obs_BMRB_diff_sub_energy
    else:
        energy += energy_if_false

    CB_obs_BMRB_diff_sub_energy = CB_obs_BMRB_diff(i, assignment)
    if CB_obs_BMRB_diff_sub_energy:
        energy += CB_obs_BMRB_diff_sub_energy
    else:
        energy += energy_if_false

    NOESY_H_dist_sub_energy = NOESY_H_dist(i, assignment, index_list)
    if NOESY_H_dist_sub_energy:
        energy += NOESY_H_dist_sub_energy
    else:
        energy += energy_if_false

    return energy


def CA_CA_prime_diff(i, assignment, index_list):  # calculates energy for difference in CA and CA-1
    try:
        peak_next = index_list[i + 1][0]

        CA_prime_shift = peak_next.get_data('CAPrimeShift')
        CA_prime_signal_noise = peak_next.get_data('HNCOCASignalNoise')

        CA_shift = assignment[0].get_data('CAShift')
        CA_signal_noise = assignment[0].get_data('HNCASignalNoise')

        if CA_prime_shift is None or CA_shift is None:
            return False

        elif abs(CA_shift - CA_prime_shift) < aa_p_d_delta and \
                CA_prime_signal_noise >= 10 and CA_signal_noise >= 10:
            return 0

        else:
            w1 = float((aa_p_d_sn_factor * ((CA_signal_noise * CA_prime_signal_noise) ** .5)) - aa_p_d_sub)
            result_new = w1 * aa_p_d * (CA_shift - CA_prime_shift)**2
            return result_new
    except IndexError:
        return False


# calculates energy for difference in CB and CB-1
def CB_CB_prime_diff(i, assignment, index_list):
    try:
        peak_next = index_list[i + 1][0]

        CB_prime_shift = index_list[i + 1][0].get_data('CBPrimeShift')
        CB_prime_signal_noise = index_list[i + 1][0].get_data('HNCOCACBSignalNoise')

        CB_shift = assignment[0].get_data('CBShift')
        CB_shift_signal_noise = assignment[0].get_data('HNCACBSignalNoise')

        if CB_prime_shift is None or CB_shift is None:
            return False

        elif abs(CB_shift - CB_prime_shift) < bb_p_d_delta and \
                CB_shift_signal_noise >= 10 and CB_prime_signal_noise >= 10:
            return 0

        else:
            w1 = float((bb_p_d_sn_factor * ((CB_shift_signal_noise * CB_prime_signal_noise) ** .5)) - bb_p_d_sub)
            result_new = w1 * bb_p_d * (CB_shift - CB_prime_shift) ** 2
            return result_new
    except IndexError:
        return False
  

# calculates energy between CA observed in spectral data, and that expected from BMRB
def CA_obs_BMRB_diff(i, assignment):
    try:
        penalty = tuple()
        CA_shift = assignment[0].get_data('CAShift')
        CA_BMRB = assignment[1].get_data('CAExpected')
        CA_BMRB_SD = assignment[1].get_data('CAExpectedSD')
        CA_signal_noise = assignment[0].get_data('HNCASignalNoise')
        if CA_shift is None:
            return False
        result = abs(CA_shift - CA_BMRB)
        if result <= CA_BMRB_SD:
            penalty = 0
        elif result > CA_BMRB_SD:
            penalty = a_bmrb_sn_factor * (CA_signal_noise ** .5) * a_bmrb * (result - CA_BMRB_SD) * (result + CA_BMRB_SD)
        return penalty
    except IndexError:
        return False


# calculates energy between CB observed in spectral data, and that expected from BMRB
def CB_obs_BMRB_diff(i, assignment):
    try:
        penalty = tuple()
        CB_shift = assignment[0].get_data('CBShift')
        CB_BMRB = assignment[1].get_data('CBExpected')
        CB_BMRB_SD = assignment[1].get_data('CBExpectedSD')
        CB_signal_noise = assignment[0].get_data('HNCACBSignalNoise')
        if CB_shift is None:
            return False
        result = abs(CB_shift - CB_BMRB)
        if result <= CB_BMRB_SD:
            penalty = 0
        elif result > CB_BMRB_SD:
            penalty = b_bmrb_sn_factor * (CB_signal_noise ** .5) * b_bmrb * (result - CB_BMRB_SD) * (result + CB_BMRB_SD)
        return penalty
    except IndexError:
        return False
    
# gets nearby H values for an assignment from noesy data. Then evaluates which residue the assignment
# corresponds to and uses the structural data to find which amino acids are closeby. Then for each of
# the amino acids that are close by, it gets their H shift by looking at their assignments.

# Essentially this generates a list of H values that we expect to see and compares it to the noesy values
# that are actually seen.


def NOESY_H_dist(i, assignment, index_list):
    residue = assignment[1]
    peak = assignment[0]

    closeby_residue_list = residue.get_data('closebyAminoAcids')
    closeby_residue_list_dist = residue.get_data('closebyAminoAcidsDist')

    if len(closeby_residue_list) == 0:
        return False

    nearby_H_shift_list = peak.get_data('NearbyHShift')
    nearby_H_shift_list_signal_noise = peak.get_data('NOESYSignalNoise')

    if len(nearby_H_shift_list) == 0:
        return False

    weight = float()

    for j, closeby in enumerate(closeby_residue_list):
        closeby_dist = closeby_residue_list_dist[j]
        closeby_assignment = index_list[closeby - 1]
        closeby_assignment_peak = closeby_assignment[0]
        closeby_peak_H_shift = closeby_assignment_peak.get_data('TROSYHShift')

        if closeby_peak_H_shift is None:
            continue

        '''
        variables so far: residue, peak, closeby_residue_list(&dist) iter: , nearbyHShiftList (& SN) 
        closebyAAnumber (&index & dist), closebyPeak (&TROSY Hshift)
        '''

        NOESY_TROSY_H_diff = 5
        index = int()
        for k, nearby_H_shift in enumerate(nearby_H_shift_list):
            if abs(nearby_H_shift - closeby_peak_H_shift) < NOESY_TROSY_H_diff:
                NOESY_TROSY_H_diff = abs(nearby_H_shift - closeby_peak_H_shift)
                index = k

        nearby_H_shift_signal_noise = nearby_H_shift_list_signal_noise[index]

        dist_factor = dist_factor_forumla(closeby_dist)
        SN_factor = sn_factor_formula(nearby_H_shift_signal_noise)
        if NOESY_TROSY_H_diff < noesy_perfect_match_threshold:
            subWeight = npmt_penatly * dist_factor * SN_factor
        elif NOESY_TROSY_H_diff < noesy_semi_perfect_match_theshold:
            subWeight = nspmt_penalty * dist_factor * SN_factor
        else:
            subWeight = n_no_match_penalty * dist_factor * SN_factor

        weight += subWeight

    return weight


'''
# CaCaPrimedistance
for i in range (0,50):
     if CA_CA_prime_diff(i,index_list[i]):
         print(CA_CA_prime_diff(i,index_list[i])[0],CA_CA_prime_diff(i,index_list[i])[1], sep=',', end=',')
         print (index_list[i][0].peak_properties['peakNumber'],
                index_list[i][0].peak_properties['HNCAHShift'],
                index_list[i][0].peak_properties['HNCANShift'],
                index_list[i][0].peak_properties['CAShift'],
                index_list[i][0].peak_properties['HNCASignalNoise'],
                sep=',', end=',')
         print (index_list[i+1][0].peak_properties['peakNumber'],
                index_list[i+1][0].peak_properties['HNCOCAHShift'],
                index_list[i+1][0].peak_properties['HNCOCANShift'],
                index_list[i+1][0].peak_properties['CAPrimeShift'],
                index_list[i+1][0].peak_properties['HNCOCASignalNoise'],
                sep=',', end = '\n')

# CbCbPrimedistance
for i in range (0,9):
    print(CB_CB_prime_diff(i,index_list[i]))
    print (index_list[i][0].peak_properties,index_list[i][1].peak_properties)

# BMRB
for i in range (0,9):
    print(CA_obs_BMRB_diff(i,index_list[i]))
    print (index_list[i][0].peak_properties,index_list[i][1].peak_properties)

for i in range (0,9):
    print(CB_obs_BMRB_diff(i,index_list[i]))
    print (index_list[i][0].peak_properties,index_list[i][1].peak_properties)

# Noesy
for i in range (0,9):
    print(NOESY_H_dist(i,index_list[i]))
    print (index_list[i][0].peak_properties,index_list[i][1].peak_properties)
'''
