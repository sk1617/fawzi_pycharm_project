import pandas as pd
import plotly.express as px
df = pd.read_table('fto_assignment_files/Full_length/Run/assigned_peaks_res.out', sep=' +')
print(df.head())

residue_count, assignment_count, low_count, med_count, high_count = 0, 0, 0, 0, 0
delta_ca_list, delta_cb_list = [], []
for row_number in range(len(df)):
    reliability = df.loc[row_number, 'reli']
    if 'F' in reliability:
        df.loc[row_number, 'reli'] = df.loc[row_number, 'reli'].replace('F', '')
        reliability = df.loc[row_number, 'reli']
    if reliability == '-':
        residue_count += 1
    else:
        residue_count += 1
        assignment_count += 1
        if 'L' in reliability:
            low_count += 1
        elif 'M' in reliability:
            med_count += 1
        elif 'H' in reliability:
            high_count += 1
    if df.loc[row_number, 'CA-1'] != '-' and df.loc[row_number - 1, 'CA'] != '-':
        delta_ca_list.append(abs(float(df.loc[row_number, 'CA-1']) - float(df.loc[row_number - 1, 'CA'])))
    else:
        delta_ca_list.append(None)
    if df.loc[row_number, 'CB-1'] != '-' and df.loc[row_number - 1, 'CB'] != '-':
        delta_cb_list.append(abs(float(df.loc[row_number, 'CB-1']) - float(df.loc[row_number - 1, 'CB'])))
    else:
        delta_cb_list.append(None)
df['Delta CA'] = delta_ca_list
df['Delta CB'] = delta_cb_list

print('Out of a total {} residues, {} are assigned, of which {} are High, {} are Med, and {} are Low'.format(
    residue_count, assignment_count, high_count, med_count, low_count))
fig_delta_ca = px.histogram(df, x='Delta CA', color='reli', nbins=100, range_x=[0, 10])
fig_delta_ca.show()
fig_delta_cb = px.histogram(df, x='Delta CB', color='reli', nbins=100, range_x=[0, 10])
fig_delta_cb.show()




