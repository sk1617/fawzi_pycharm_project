from Bio.PDB import PDBParser

# create parser
parser = PDBParser()

# read structure from file. The pdb file should be in the same directory as the python file
structure = parser.get_structure('Ailin\'s protein', '2n4p.pdb')

# model is the first structure of the pdb
model = structure[0]
# for Ailin's protein, there is only one chain.
chain = model['A']

# chains are made of class Residue, residues are made of class Atom

c = 0
for residue1 in chain:
    if 13 <= residue1.full_id[3][1]:
        for residue2 in chain:
            if 13 <= residue2.full_id[3][1]:
                if residue1 != residue2:
                    # compute distance between CA atoms
                    try:
                        # the subtract func has been superseded to find distances between two atoms
                        distance = residue1['H'] - residue2['H']
                    except KeyError:
                        # no H atom for Pro
                        continue
                    if distance < 15:  # change this number for smaller/larger dist
                        print(residue1.get_resname(), residue1.full_id[3][1] - 12, 'HN',
                              residue2.get_resname(), residue2.full_id[3][1] - 12, 'HN',
                              distance, sep='   ')
