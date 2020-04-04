#%%
from EnergyFunctions import *
from DataImporting import *
import plotly.express as px
import plotly.graph_objects as go
import BMRB

#%%
# Old CaCa'
CA_prime_signal_noise = 10
CA_signal_noise = 10
def ca(delta):
    if abs(delta) < aa_p_d_delta and \
            CA_prime_signal_noise >= 10 and CA_signal_noise >= 10:
        return 0.1

    else:
        w1 = float((aa_p_d_sn_factor * ((CA_signal_noise * CA_prime_signal_noise) ** .5)) - aa_p_d_sub)
        result_new = w1 * aa_p_d * delta**2
        return result_new
x = [i/100 for i in range(400)]
old_ca_delta_plt = px.line(x=x, y=[ca(i) for i in x], title='OLD Ca-CaPrime Energy by Delta ppm')
old_ca_delta_plt.show()

#%%
# New CaCa'
aa_p_d = 1875
def ca_new(delta):
    weight = aa_p_d * delta**2 + 0.000001
    return weight
x = [i/100 for i in range(100)]
new_ca_delta_plt = px.line(x=x, y=[ca_new(i) for i in x], title='NEW Ca-CaPrime Energy by Delta ppm')
new_ca_delta_plt.show()

#%%
CB_prime_signal_noise = 10
CB_shift_signal_noise = 10
def cb(delta):
    if abs(delta) < bb_p_d_delta and \
            CB_shift_signal_noise >= 10 and CB_prime_signal_noise >= 10:
        return 0.1
    else:
        w1 = float((bb_p_d_sn_factor * ((CB_shift_signal_noise * CB_prime_signal_noise) ** .5)) - bb_p_d_sub)
        result_new = w1 * bb_p_d * delta**2
        return result_new
old_cb_delta_plt = px.line(x=x, y=[cb(i) for i in x], title='OLD CbCbPrime energy by Delta ppm')
old_cb_delta_plt.show()

#%%
bb_p_d = 1875
def cb_new(delta):
    weight = bb_p_d * delta**2 + 1e-6
    return weight
new_cb_delta_plt = px.line(x=x, y=[cb_new(i) for i in x], title='NEW CbCbPrime Energy by Delta ppm')
new_cb_delta_plt.show()

#%%

# def BMRB_diff(i, peak_assignment_index, peak_list, residue_list, Delta_List, append=False, use_prime_data=False):
def BMRB_diff_new(ca_delta, cb_delta):
    sub_energy = float()
    CA_BMRB_SD = 1
    CB_BMRB_SD = 1
    ca_energy = 158 * ca_delta**2 if ca_delta else energy_if_false
    cb_energy = 158 * cb_delta**2 if cb_delta else energy_if_false
    sub_energy = ca_energy + cb_energy
    return sub_energy

x, y, z = [], [], []
a = [i/10 for i in range(30)]
b = [i/10 for i in range(30)]
for i in a:
    for j in b:
        x.append(i)
        y.append(j)
        z.append(BMRB_diff_new(i, j))
new_bmrb_fig = px.scatter_3d(x=x, y=y, z=z, title='NEW BMRB')
new_bmrb_fig.show()

#%%

# def BMRB_diff(i, peak_assignment_index, peak_list, residue_list, Delta_List, append=False, use_prime_data=False):
def BMRB_diff_old(ca_delta, cb_delta):
    sub_energy = float()
    CA_BMRB_SD = 1
    CB_BMRB_SD = 1
    if ca_delta > CA_BMRB_SD and cb_delta > CB_BMRB_SD:
        sub_energy = ((ca_delta - CA_BMRB_SD) * (ca_delta + CA_BMRB_SD) *
                      (cb_delta - CB_BMRB_SD) * (cb_delta + CB_BMRB_SD)) * bmrb
    elif ca_delta < CA_BMRB_SD and cb_delta > CB_BMRB_SD:
        sub_energy = ((cb_delta - CB_BMRB_SD) * (cb_delta + CB_BMRB_SD)) ** 2 * bmrb
    elif ca_delta > CA_BMRB_SD and cb_delta < CB_BMRB_SD:
        sub_energy = ((ca_delta - CA_BMRB_SD) * (ca_delta + CA_BMRB_SD)) ** 2 * bmrb
    elif ca_delta < CA_BMRB_SD and cb_delta < CB_BMRB_SD:
        sub_energy = 0.1

    return sub_energy

x, y, z = [], [], []
a = [i/10 for i in range(30)]
b = [i/10 for i in range(30)]
for i in a:
    for j in b:
        x.append(i)
        y.append(j)
        z.append(BMRB_diff_old(i, j))
old_bmrb_fig = px.scatter_3d(x=x, y=y, z=z, title='OLD BMRB')
old_bmrb_fig.show()