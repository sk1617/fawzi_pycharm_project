import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
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


from DataProcessing import main_data_processing, distance_formula
outside_assignment_df = pd.read_table('fto_assignment_files/Full_length/Run/assigned_peaks_res.out', sep=' +')
outside_peaks_df = pd.read_table('fto_assignment_files/Full_length/CS17_FTO', sep=' +', index_col=False)
peak_list, residue_list = main_data_processing()

print('NH')
total_rows = len(outside_peaks_df)
print(total_rows)
print('CA')
print(total_rows - sum(outside_peaks_df['CA'].str.count('-')))
print('CB')
print(total_rows - sum(outside_peaks_df['CB'].str.count('-')))
print('CA-1')
print(total_rows - sum(outside_peaks_df['CA-1'].str.count('-')))
print('CB-1')
print(total_rows - sum(outside_peaks_df['CB-1'].str.count('-')))
print('CO-1')
print(total_rows - sum(outside_peaks_df['CO-1'].str.count('-')))

h_old = []
n_old = []
h_new = [float(x) for x in outside_peaks_df['H']]
n_new = [float(x) for x in outside_peaks_df['N']]

for peak in peak_list:
    h_shift = peak.get_data('TROSYHShift')
    n_shift = peak.get_data('TROSYNShift')
    if h_shift == None or n_shift is None:
        continue
    else:
        h_old.append(float(h_shift))
        n_old.append(float(n_shift))

fig = go.Figure()
fig.add_trace(go.Scatter(x=h_new, y=n_new, mode='markers', name='new'))
fig.add_trace(go.Scatter(x=h_old, y=n_old, mode='markers', name='old'))
fig.update_layout(title='Comparing new and old NH graphs', xaxis_title='H (ppm)', yaxis_title='N (ppm)')
fig.show()

import chart_studio.plotly as py
import chart_studio
chart_studio.tools.set_credentials_file(username='sk1617', api_key='htVpHm0fYKP0sT2H5vHu')
#
py.iplot(fig)