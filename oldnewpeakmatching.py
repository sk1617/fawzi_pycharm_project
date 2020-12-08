import pandas as pd
from DataProcessing import main_data_processing, distance_formula
import plotly.graph_objects as go
outside_assignment_df = pd.read_table('fto_assignment_files/Full_length/Run/assigned_peaks_res.out', sep=' +')
outside_peaks_df = pd.read_table('fto_assignment_files/Full_length/CS17_FTO', sep=' +', index_col=False)
peak_list, residue_list = main_data_processing()

outside_peaks_df['N_'] = [float(x)+.23 for x in outside_peaks_df['N']]
outside_peaks_df['H_'] = [float(x)+.08 for x in outside_peaks_df['H']]


# 3D HN (CA | CB | CA' | CB')
h_old, n_old = [], []
ca_old, ca_prime_old = [], []
cb_old, cb_prime_old = [], []
h_new = outside_peaks_df['H_']
n_new = outside_peaks_df['N_']
ca_new = [float(x) if x != '-' else 0 for x in outside_peaks_df['CA']]
ca_prime_new = [float(x) if x != '-' else 0 for x in outside_peaks_df['CA-1']]
cb_new = [float(x) if x!= '-' else 0 for x in outside_peaks_df['CB']]
cb_prime_new = [float(x) if x != '-' else 0 for x in outside_peaks_df['CB-1']]

for peak in peak_list:
    h_shift = peak.get_data('TROSYHShift')
    n_shift = peak.get_data('TROSYNShift')
    ca_shift = peak.get_data('CAShift')
    ca_prime_shift = peak.get_data('CAPrimeShift')
    cb_shift = peak.get_data('CBShift')
    cb_prime_shift = peak.get_data('CBPrimeShift')
    if h_shift == None or n_shift is None:
        continue
    else:
        h_old.append(float(h_shift))
        n_old.append(float(n_shift))
        ca_old.append(ca_shift) if ca_shift is not None else ca_old.append(0)
        ca_prime_old.append(ca_prime_shift) if ca_prime_shift is not None else ca_prime_old.append(0)
        cb_old.append(cb_shift) if cb_shift is not None else cb_old.append(0)
        cb_prime_old.append(cb_prime_shift) if cb_prime_shift is not None else cb_prime_old.append(0)


fig = go.Figure()
# fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=ca_prime_old, mode='markers', name='old'))
# fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=ca_prime_new, mode='markers', name='new'))
# fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=ca_old, mode='markers', name='old'))
# fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=ca_new, mode='markers', name='new'))
# fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=cb_old, mode='markers', name='old'))
# fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=cb_new, mode='markers', name='new'))
fig.add_trace(go.Scatter3d(x=h_old, y=n_old, z=cb_prime_old, mode='markers', name='old'))
fig.add_trace(go.Scatter3d(x=h_new, y=n_new, z=cb_prime_new, mode='markers', name='new'))

fig.update_layout(title='Comparing new and old NH-CB Prime graphs', xaxis_title='H (ppm)', yaxis_title='N (ppm)')
fig.show()



# 2D old and new NH peaks
fig = go.Figure()
fig.add_trace(go.Scatter(x=h_old, y=n_old, mode='markers', name='old'))
fig.add_trace(go.Scatter(x=h_new, y=n_new, mode='markers', name='new'))
fig.update_layout(title='Comparing new and old NH graphs', xaxis_title='H (ppm)', yaxis_title='N (ppm)')
fig.show()


# Graphing the spread of CB/A and CB/A Prime shifts of the new and old data as a histogram
cb_old, cb_prime_old = [], []
ca_old, ca_prime_old = [], []
co_old = []
for peak in peak_list:
    cb, cb_prime = peak.get_data('CBShift'), peak.get_data('CBPrimeShift')
    ca, ca_prime = peak.get_data('CAShift'), peak.get_data('CAPrimeShift')
    co_prime = peak.get_data('COPrimeShift')
    if cb is not None:
        cb_old.append(cb)
    if cb_prime is not None:
        cb_prime_old.append(cb_prime)
    if ca is not None:
        ca_old.append(ca)
    if ca_prime is not None:
        ca_prime_old.append(ca_prime)
    if co_prime is not None:
        co_old.append(co_prime)

cb_new = [float(x) if x!= '-' else '-' for x in outside_peaks_df['CB']]
cb_prime_new = [float(x) if x != '-' else '-' for x in outside_peaks_df['CB-1']]
ca_new = [float(x) if x!= '-' else '-' for x in outside_peaks_df['CA']]
ca_prime_new = [float(x) if x != '-' else '-' for x in outside_peaks_df['CA-1']]
co_new = [float(x) if x != '-' else '-' for x in outside_peaks_df['CO-1']]


# plotting Histograms
if True:
    figCB = go.Figure()
    figCB.add_trace(go.Histogram(x = cb_old, name='old'))
    figCB.add_trace(go.Histogram(x = cb_new, name='new'))
    figCB.update_layout(barmode = 'overlay', title_text = 'New CB and Old CB Shifts Histograms Overlaid',
                        xaxis_title_text='CB Shift (ppm)', yaxis_title_text='count')
    figCB.update_traces(opacity=0.75)
    figCB.show()

    figCBPrime = go.Figure()
    figCBPrime.add_trace(go.Histogram(x = cb_prime_old, name='old'))
    figCBPrime.add_trace(go.Histogram(x = cb_prime_new, name='new'))
    figCBPrime.update_layout(barmode = 'overlay', title_text = 'New CB Prime and Old CB Prime Shifts Histograms Overlaid',
                             xaxis_title_text='CB Prime Shift (ppm)', yaxis_title_text='count')
    figCBPrime.update_traces(opacity=0.75)
    figCBPrime.show()

    figCA = go.Figure()
    figCA.add_trace(go.Histogram(x = ca_old, name='old'))
    figCA.add_trace(go.Histogram(x = ca_new, name='new'))
    figCA.update_layout(barmode = 'overlay', title_text = 'New CA and Old CA Shifts Histograms Overlaid',
                        xaxis_title_text='CA Shift (ppm)', yaxis_title_text='count')
    figCA.update_traces(opacity=0.75)
    figCA.show()

    figCAPrime = go.Figure()
    figCAPrime.add_trace(go.Histogram(x = ca_prime_old, name='old'))
    figCAPrime.add_trace(go.Histogram(x = ca_prime_new, name='new'))
    figCAPrime.update_layout(barmode = 'overlay', title_text = 'New CA Prime and Old CA Prime Shifts Histograms Overlaid',
                             xaxis_title_text='CA Prime Shift (ppm)', yaxis_title_text='count')
    figCAPrime.update_traces(opacity=0.75)
    figCAPrime.show()

    figCOPrime = go.Figure()
    figCOPrime.add_trace(go.Histogram(x = co_old, name='old'))
    figCOPrime.add_trace(go.Histogram(x = co_new, name='new'))
    figCOPrime.update_layout(barmode = 'overlay', title_text = 'New CO Prime and Old CO Prime Shifts Histograms Overlaid',
                             xaxis_title_text='CO Prime Shift (ppm)', yaxis_title_text='count')
    figCOPrime.update_traces(opacity=0.75)
    figCOPrime.show()

# 2D N- (CA, CA', CB, CB') graphs
figNCA = go.Figure()
figNCA.add_trace(go.Scatter(x=n_old, y=ca_old, mode='markers', name='old'))
figNCA.add_trace(go.Scatter(x=n_new, y=ca_new, mode='markers', name='new'))
figNCA.update_layout(title='Comparing new and old CA - N graphs', xaxis_title='N (ppm)', yaxis_title='CA (ppm)')
figNCA.show()

figNCAPrime = go.Figure()
figNCAPrime.add_trace(go.Scatter(x=n_old, y=ca_prime_old, mode='markers', name='old'))
figNCAPrime.add_trace(go.Scatter(x=n_new, y=ca_prime_new, mode='markers', name='new'))
figNCAPrime.update_layout(title='Comparing new and old CA Prime - N graphs', xaxis_title='N (ppm)', yaxis_title='CA Prime (ppm)')
figNCAPrime.show()

figNCB = go.Figure()
figNCB.add_trace(go.Scatter(x=n_old, y=cb_old, mode='markers', name='old'))
figNCB.add_trace(go.Scatter(x=n_new, y=cb_new, mode='markers', name='new'))
figNCB.update_layout(title='Comparing new and old CB - N graphs', xaxis_title='N (ppm)', yaxis_title='CB (ppm)')
figNCB.show()

figNCBPrime = go.Figure()
figNCBPrime.add_trace(go.Scatter(x=n_old, y=cb_prime_old, mode='markers', name='old'))
figNCBPrime.add_trace(go.Scatter(x=n_new, y=cb_prime_new, mode='markers', name='new'))
figNCBPrime.update_layout(title='Comparing new and old CB Prime- N graphs', xaxis_title='N (ppm)', yaxis_title='CB Prime (ppm)')
figNCBPrime.show()

