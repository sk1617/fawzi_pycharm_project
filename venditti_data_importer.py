from DataImporting import sequence_list
fh = open('venditti_assignment.txt')
import re

file = fh.read()
ls = re.findall(r"label \d*N-H\nxy \d*.\d*,\d*.\d*", file)

venditti_assignment = list(range(0, len(sequence_list)))

for i in ls:
    index = int(re.search(r"\d*N-H", i)[0][:-3])
    if index > 1000:
        break
    h_shift = re.search(r",\d\.\d*", i)[0][1:]
    n_shift = re.search(r" \d\d\d\.\d*", i)[0][1:]
    assignment_tuple = (index, h_shift, n_shift)
    venditti_assignment[index - 1] = assignment_tuple