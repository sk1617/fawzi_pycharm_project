from DataImporting import *

# @profile
def main_energy_function(a, b, old_index_list, new_index_list, peak_list, residue_list, energy_if_false, Delta_List):
    old_energy, new_energy = 0.0, 0.0
    energy_list = [old_energy, new_energy]
    for i, index_list in enumerate([old_index_list, new_index_list]):
        list_for_prime = [min(a) - 1] + [x for x in a] + [min(b) - 1] + [x for x in b]  # move this out of the loop
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

        CA_shift = peak.get_data('CAShift')
        CA_prime_shift = peak_next.get_data('CAPrimeShift')

        if CA_prime_shift is None or CA_shift is None:
            return False
        else:
            result_new = aa_p_d * (CA_shift - CA_prime_shift)**2 + 1e-6
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
        CB_shift = peak.get_data('CBShift')

        if CB_prime_shift is None or CB_shift is None:
            return False
        else:
            result_new = bb_p_d * (CB_shift - CB_prime_shift) ** 2 + 1e-6
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
    if residue.get_data('residue_name') == 'G' and CB_shift is not None:
        return 10000
    CB_BMRB = residue.get_data('CBExpected')
    CB_BMRB_SD = residue.get_data('CBExpectedSD')

    ca_delta = (CA_shift - CA_BMRB)/CA_BMRB_SD if CA_shift else False
    try:
        cb_delta = (CB_shift - CB_BMRB)/CB_BMRB_SD if CB_shift else False
    except ZeroDivisionError:
        print(CB_shift, CB_BMRB, CA_BMRB_SD, residue.residue_properties, peak.peak_properties)
        raise ZeroDivisionError

    ca_energy = bmrb_ca * ca_delta**2 if ca_delta else energy_if_false
    cb_energy = bmrb_cb * cb_delta**2 if cb_delta else energy_if_false
    sub_energy = ca_energy + cb_energy + 1e-6

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

    peak = peak_list[pai]  # gets peak
    residue = residue_list[i]  # gets residue
    assert i + 1 == residue.get_data('index')

    closeby_residue_list = residue.get_data('closebyAminoAcids')  # gets list of closeby residues from residue
    nearby_H_shift_list = peak.get_data('NearbyHShift')  # gets list of nearby H shifts from peak
    if len(closeby_residue_list) == 0:  # used
        return False
    if len(nearby_H_shift_list) == 0:  # used
        return False

    closeby_residues_peak_list = [peak_list[index_list[closeby - 1]] for closeby in closeby_residue_list]
    closeby_residues_peak_h_shift = [peak.get_data('TROSYHShift') for peak in closeby_residues_peak_list if
                                     peak.get_data('TROSYHShift') is not None]
    if len(closeby_residues_peak_h_shift) == 0:
        return False

    weight = 1e-6
    for noesy_shift in nearby_H_shift_list:
        closest_match_delta = min([abs(i - noesy_shift) for i in closeby_residues_peak_h_shift])
        weight += noesy_weight * closest_match_delta**2

    if should_append_DL or append:
        Delta_List.append(('NOESY', weight, [], []))
    return weight


