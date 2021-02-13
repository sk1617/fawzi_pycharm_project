
import plotly.express as px
import plotly.graph_objects as go
from MMAl import *
import MMAl
import pandas as pd

peak_list, residue_list = main_data_processing()

#%%

new_noesy_df = pd.read_csv('/Users/sohamkale/Downloads/3D_NOESY_for_Soham-selected/NOESY_HN_CS_adjusted_to_assn', sep=' +')
new_noesy_df

#%%

h_old_ls, n_old_ls = [], []
for peak in peak_list:
    h_old = peak.get_data('TROSYHShift')
    n_old = peak.get_data('TROSYNShift')
    if h_old is not None:
        h_old_ls.append(h_old)
        n_old_ls.append(n_old)

#%%

h_new_ls = list(new_noesy_df['w3(HN)']+.08)
n_new_ls = list(new_noesy_df['w2(N)']+.23)

#%%

# add figure
fig1 = go.Figure()
fig1.add_trace(go.Scatter(x=h_old_ls, y=n_old_ls, mode='markers', name='Early 2020 TROSY (N-H)', marker_color='green', opacity=0.75))
fig1.add_trace(go.Scatter(x=h_new_ls, y=n_new_ls, mode='markers', name='2021 NOESY (N-NH)', marker_color='red', opacity=.75))
fig1.update_yaxes(autorange="reversed")
fig1.update_xaxes(autorange="reversed")
fig1.update_layout(yaxis_title_text='Delta N ppm',
                   xaxis_title_text='Delta H ppm', title_text = 'New NOESY Matching Ability')
fig1.show()