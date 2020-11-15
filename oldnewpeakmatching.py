import pandas as pd
from DataProcessing import main_data_processing, distance_formula
import plotly.graph_objects as go
outside_assignment_df = pd.read_table('fto_assignment_files/Full_length/Run/assigned_peaks_res.out', sep=' +')
outside_peaks_df = pd.read_table('fto_assignment_files/Full_length/CS17_FTO', sep=' +', index_col=False)
peak_list, residue_list = main_data_processing()

outside_peaks_df['N_'] = [float(x)+.23 for x in outside_peaks_df['N']]
outside_peaks_df['H_'] = [float(x)+.08 for x in outside_peaks_df['H']]

h_old, n_old = [], []
ca_old, ca_prime_old = [], []
h_new = outside_peaks_df['H_']
n_new = outside_peaks_df['N_']
ca_new = [float(x) if x != '-' else 0 for x in outside_peaks_df['CA']]
ca_prime_new = [float(x) if x != '-' else 0 for x in outside_peaks_df['CA-1']]

for peak in peak_list:
    h_shift = peak.get_data('TROSYHShift')
    n_shift = peak.get_data('TROSYNShift')
    ca_shift = peak.get_data('CAShift')
    ca_prime_shift = peak.get_data('CAPrimeShift')
    if h_shift == None or n_shift is None:
        continue
    else:
        h_old.append(float(h_shift))
        n_old.append(float(n_shift))
        ca_old.append(ca_shift) if ca_shift is not None else ca_old.append(0)
        ca_prime_old.append(ca_prime_shift) if ca_prime_shift is not None else ca_prime_old.append(0)

fig = go.Figure()
fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=ca_prime_old, mode='markers', name='old'))
# fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=ca_old, mode='markers', name='old'))
fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=ca_prime_new, mode='markers', name='new'))
# fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=ca_new, mode='markers', name='new'))

fig.update_layout(title='Comparing new and old NH graphs', xaxis_title='H (ppm)', yaxis_title='N (ppm)')
fig.show()

fig2d = go.Figure
fig.add_trace



