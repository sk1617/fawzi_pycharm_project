from tpe_DataProcessing import *
from EnergyFunctions import *
peak_list, residue_list = main_data_processing()
import pandas as pd

perfect_list = [x for x in range(74)]

it = sequence_list.__iter__()
w = 80
v = 0
for u, thing in enumerate(it):
    if thing == 'P':
        v += 1
        perfect_list.insert(u, w - v)

print(perfect_list)
delta_list = eval_energy(perfect_list, peak_list, residue_list, 300.)[1]

# delta_list_iter, perfect_list, sequence_list, peak list, residue list

int_iter = iter([x for x in range(500)])
delta_list_iter = iter(delta_list)
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
        for i in closebyAA:
            CBAAHShift.append(peak_list[perfect_list[i - 1]].get_data('TROSYHShift'))
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

print(data_frame)

writer = pd.ExcelWriter('output for lower NRG.xlsx')
data_frame.to_excel(writer,'low NRG list 2')
