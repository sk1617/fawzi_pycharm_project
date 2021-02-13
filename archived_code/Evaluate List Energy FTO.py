from DataProcessing import *
from archived_code.Pattern_Analysis import mode_list
from EnergyFunctions import *
import pandas as pd
peak_list, residue_list = main_data_processing()


# the mode has things like None and NA. must be sanitized before putting into eval_energy of DP
mode_list_sanitized = list()
for assignment in mode_list:
    if type(assignment) == int:
        mode_list_sanitized.append(assignment)
    else:
        mode_list_sanitized.append(len(mode_list) - 1)

# evaluate the energies
energy, info_list = eval_energy(mode_list_sanitized, peak_list, residue_list, energy_if_false)

# putting it into Data Frame
int_iter = iter([x for x in range(500)])
delta_list_iter = iter(info_list)
line = next(delta_list_iter)
list = []
while True:
    try:
        if type(line[0]) is int:
            assert next(int_iter) == line[0]
            pai = line[1]
            index = line[0]
        # residue
        residue_name = sequence_list[index]
        residue = residue_list[index]
        CAExp = residue.get_data('CAExpected')
        CBExp = residue.get_data('CBExpected')
        closebyAA = residue.get_data('closebyAminoAcids')
        CBAAHShift = []
        peak = peak_list[pai]
        CAShift = peak.get_data('CAShift')
        CBShift = peak.get_data('CBShift')
        CAPrimeShift = peak.get_data('CAPrimeShift')
        CBPrimeShift = peak.get_data('CBPrimeShift')
        nearbyHShift = peak.get_data('NearbyHShift')

        line = next(delta_list_iter)
        if line[0] == 'CACAPrime':
            CANRG = line[1]
            CADel = line[2]
            line = next(delta_list_iter)
        else:
            CANRG = 200
            CADel = 0

        if line[0] == 'CBCBPrime':
            CBNRG = line[1]
            CBDel = line[2]
            line = next(delta_list_iter)
        else:
            CBNRG = 200
            CBDel = 0

        if line[0] == 'CACBExp':
            CACBExpNRG = line[1]
            CAExpDel = line[2]
            CBExpDel = line[3]
            line = next(delta_list_iter)
        else:
            CACBExpNRG = 200
            CAExpDel = 0
            CBExpDel = 0

        if line[0] == 'CACBExpPrime':
            CACBExpPrimeNRG = line[1]
            CAExpPrimeDel = line[2]
            CBExpPrimeDel = line[3]
            line = next(delta_list_iter)
        else:
            CACBExpPrimeNRG = 200
            CAExpPrimeDel = 0
            CBExpPrimeDel = 0

        if line[0] == 'NOESY':
            NOESYNRG = line[1]
            line = next(delta_list_iter)
        else:
            NOESYNRG = 200

        list.append([pai, residue_name, CAShift, CAPrimeShift, CADel, CANRG, CBShift,
                           CBPrimeShift, CBDel, CBNRG,
                           CAExp, CAExpDel, CBExp, CBExpDel, CACBExpNRG, CACBExpPrimeNRG,
                           NOESYNRG, closebyAA, CBAAHShift, nearbyHShift])

    except StopIteration:
        break

data_frame = pd.DataFrame(list, columns=['PAI', 'Residue',
                                  'CAShift', 'CAPrimeShift', 'CADel', 'CANRG', 'CBShift', 'CBPrimeShift', 'CBDel', 'CBNRG',
                                  'CAExp', 'CAExpDel', 'CBExp', 'CBExpDel', 'CACBExpNRG', 'CACBExpPrimeNRG',
                                  'NOESYNRG', 'closebyAA', 'CBAA\'s Hshift', 'nearbyHShift'])

writer = pd.ExcelWriter('1e8_25_list_eval.xlsx')
data_frame.to_excel(writer,'SHEEET')