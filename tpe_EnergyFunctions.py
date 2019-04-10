from DataImporting import *

# @profile
def main_energy_function(a, b, old_index_list, new_index_list, peak_list, residue_list, energy_if_false, Delta_List):
    old_energy, new_energy = 0.0, 0.0
    energy_list = [old_energy, new_energy]
    for i, index_list in enumerate([old_index_list, new_index_list]):
        list_for_prime = [min(a) - 1] + [x for x in a] + [min(b) - 1] + [x for x in b] # move this out of the loop
        list_for_expected = [x for x in a] + [x for x in b]

        # CA/CA Prime
        for j in (list_for_prime):
            if j < 0: continue
            en = CA_CA_prime_diff(j, index_list[j], index_list, peak_list, Delta_List)
            energy_list[i] += en if en else energy_if_false

        # CB/CB Prime
        for j in (list_for_prime):
            if j < 0: continue
            en = CB_CB_prime_diff(j, index_list[j], index_list, peak_list, Delta_List)
            energy_list[i] += en if en else energy_if_false

        # expected for CA and CB
        for j in (list_for_expected):
            energy_list[i] += BMRB_diff(j, index_list[j], peak_list, residue_list, Delta_List)

        for j in [min(a)] + [min(b)]:
            energy_list[i] += BMRB_diff(j, index_list[j], peak_list, residue_list, Delta_List, use_prime_data=True)

        # noesy data
        for j in (list_for_expected):
            energy_list[i] += NOESY_H_dist(j, index_list[j], index_list, peak_list, residue_list, Delta_List)

            residue = residue_list[j]
            closeby_residue_list = residue.get_data('closebyAminoAcids')

            for k in closeby_residue_list:
                energy_list[i] += NOESY_H_dist(k - 1, index_list[k - 1], index_list, peak_list, residue_list, Delta_List)

    delta = energy_list[1] - energy_list[0]
    return delta


def eval_energy(index_list, peak_list, residue_list, energy_if_false):
    info_list = []
    energy = float()
    for i, pai in enumerate(index_list):
        info_list.append((i, pai))
        en = CA_CA_prime_diff(i, pai, index_list, peak_list, info_list, append=True)
        energy += en if en else energy_if_false
        en = CB_CB_prime_diff(i, pai, index_list, peak_list, info_list, append=True)
        energy += en if en else energy_if_false
        en = BMRB_diff(i, pai, peak_list, residue_list, info_list, append=True)
        energy += en if en else energy_if_false
        en = BMRB_diff(i, pai, peak_list, residue_list, info_list, append=True, use_prime_data=True)
        energy += en if en else energy_if_false
        en = NOESY_H_dist(i, pai, index_list, peak_list, residue_list, info_list, append=True)
        energy += en if en else energy_if_false
    info_list.append((energy, index_list))
    return energy, info_list


# calculates energy for difference in CA and CA-1
# @profile
def CA_CA_prime_diff(i, peak_assignment_index, index_list, peak_list, Delta_List, append=False):
    try:
        peak = peak_list[peak_assignment_index]
        peak_next = peak_list[index_list[i + 1]]

        CA_prime_shift = peak_next.get_data('CAPrimeShift')
        CA_prime_signal_noise = peak_next.get_data('HNCOCASignalNoise')

        CA_shift = peak.get_data('CAShift')
        CA_signal_noise = peak.get_data('HNCASignalNoise')

        if CA_prime_shift is None or CA_shift is None:
            return False
        elif abs(CA_shift - CA_prime_shift) < aa_p_d_delta and \
                CA_prime_signal_noise >= 10 and CA_signal_noise >= 10:
            if should_append_DL or append:
                Delta_List.append(('CACAPrime', 0.1, CA_shift - CA_prime_shift))
            return 0.1

        else:
            w1 = float((aa_p_d_sn_factor * ((CA_signal_noise * CA_prime_signal_noise) ** .5)) - aa_p_d_sub)
            result_new = w1 * aa_p_d * (CA_shift - CA_prime_shift)**2
            if should_append_DL or append:
                Delta_List.append(('CACAPrime', result_new, CA_shift - CA_prime_shift))
            return result_new
    except IndexError:
        return False


# calculates energy for difference in CB and CB-1
def CB_CB_prime_diff(i, pk_ass_ind, index_list, peak_list, Delta_List, append=False):
    try:
        peak = peak_list[pk_ass_ind]
        peak_next = peak_list[index_list[i + 1]]

        CB_prime_shift = peak_next.get_data('CBPrimeShift')
        CB_prime_signal_noise = peak_next.get_data('HNCOCACBSignalNoise')

        CB_shift = peak.get_data('CBShift')
        CB_shift_signal_noise = peak.get_data('HNCACBSignalNoise')

        if CB_prime_shift is None or CB_shift is None:
            return False

        elif abs(CB_shift - CB_prime_shift) < bb_p_d_delta and \
                CB_shift_signal_noise >= 10 and CB_prime_signal_noise >= 10:
            if should_append_DL or append:
                Delta_List.append(('CBCBPrime', 0.1, CB_shift - CB_prime_shift))
            return 0.1

        else:
            w1 = float((bb_p_d_sn_factor * ((CB_shift_signal_noise * CB_prime_signal_noise) ** .5)) - bb_p_d_sub)
            result_new = w1 * bb_p_d * (CB_shift - CB_prime_shift) ** 2
            if should_append_DL or append:
                Delta_List.append(('CBCBPrime', result_new, CB_shift - CB_prime_shift))
            return result_new

    except IndexError:
        return False
  

# calculates energy between CA observed in spectral data, and that expected from BMRB
# @profile
def BMRB_diff(i, peak_assignment_index, peak_list, residue_list, Delta_List, append=False, use_prime_data=False):
    if i == 0 and use_prime_data:
        return energy_if_false
    sub_energy = float()
    peak = peak_list[peak_assignment_index]
    residue = residue_list[i] if not use_prime_data else residue_list[i - 1]

    CA_shift = peak.get_data('CAShift') if not use_prime_data else peak.get_data('CAPrimeShift')
    CA_BMRB = residue.get_data('CAExpected')
    CA_BMRB_SD = residue.get_data('CAExpectedSD')

    CB_shift = peak.get_data('CBShift') if not use_prime_data else peak.get_data('CBPrimeShift')
    CB_BMRB = residue.get_data('CBExpected')
    CB_BMRB_SD = residue.get_data('CBExpectedSD')

    if CA_shift is None and CB_shift is None:
        sub_energy += energy_if_false * 2
        ca_delta = 0
        cb_delta = 0

    if CA_shift is None:
        sub_energy += energy_if_false
        ca_delta = 0
    else:
        ca_delta = abs(CA_shift - CA_BMRB)

    if CB_shift is None:
        sub_energy += energy_if_false
        cb_delta = 0
    else:
        cb_delta = abs(CB_shift - CB_BMRB)

    if ca_delta > CA_BMRB_SD and cb_delta > CB_BMRB_SD:
        sub_energy = ((ca_delta - CA_BMRB_SD) * (ca_delta + CA_BMRB_SD) * \
                      (cb_delta - CB_BMRB_SD) * (cb_delta + CB_BMRB_SD)) * \
                     bmrb

    elif ca_delta < CA_BMRB_SD and cb_delta > CB_BMRB_SD:
        sub_energy = ((cb_delta - CB_BMRB_SD) * (cb_delta + CB_BMRB_SD)) ** 2 * \
                     bmrb
    elif ca_delta > CA_BMRB_SD and cb_delta < CB_BMRB_SD:
        sub_energy = ((ca_delta - CA_BMRB_SD) * (ca_delta + CA_BMRB_SD)) ** 2 * \
                     bmrb
    elif ca_delta < CA_BMRB_SD and cb_delta < CB_BMRB_SD:
        sub_energy = 0.1

    if should_append_DL or append:
        if not use_prime_data:
            Delta_List.append(('CACBExp', sub_energy, ca_delta, cb_delta))
        else:
            Delta_List.append(('CACBExpPrime', sub_energy, ca_delta, cb_delta))

    return sub_energy


# gets nearby H values for an assignment from noesy data. Then evaluates which residue the assignment
# corresponds to and uses the structural data to find which amino acids are closeby. Then for each of
# the amino acids that are close by, it gets their H shift by looking at their assignments.

# Essentially this generates a list of H values that we expect to see and compares it to the noesy values
# that are actually seen.


# @profile
def NOESY_H_dist(i, pai, index_list, peak_list, residue_list, Delta_List, append=False):

    peak = peak_list[pai]
    residue = residue_list[i]

    closeby_residue_list = residue.get_data('closebyAminoAcids')
    closeby_residue_list_dist = residue.get_data('closebyAminoAcidsDist')

    if len(closeby_residue_list) == 0:
        return False

    nearby_H_shift_list = peak.get_data('NearbyHShift')
    nearby_H_shift_list_signal_noise = peak.get_data('NOESYSignalNoise')

    if len(nearby_H_shift_list) == 0:
        return False

    weight = float()
    noesy_trosy_diff_list = []
    subweight_list =[]

    for j, closeby in enumerate(closeby_residue_list):
        closeby_dist = closeby_residue_list_dist[j]
        closeby_assignment_peak = peak_list[index_list[closeby - 1]]

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
        elif NOESY_TROSY_H_diff < noesy_semi_perfect_match_threshold:
            subWeight = nspmt_penalty * dist_factor * SN_factor
        else:
            subWeight = n_no_match_penalty * dist_factor * SN_factor

        noesy_trosy_diff_list.append(NOESY_TROSY_H_diff)
        subweight_list.append(subWeight)
        weight += subWeight

    if should_append_DL or append:
        Delta_List.append(('NOESY', weight, subweight_list, noesy_trosy_diff_list))
    return weight


