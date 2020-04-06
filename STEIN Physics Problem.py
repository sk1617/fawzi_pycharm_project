import numpy as np
import plotly.graph_objects as go
def prob(i, ep = 1.6567788e-21, s = 5.17743375e-24, kb = 1.380649e-23, N = 10):
    return np.exp((N/kb) * ((ep/i) - s))
x = [i for i in range(50, 700)]

# Changing N
plt = go.Figure()
plt.add_trace(go.Scatter(x=x, y=[prob(i, N = 10) for i in x], name = 'N = 10'))
plt.add_trace(go.Scatter(x=x, y=[prob(i, N = 20) for i in x], name = 'N = 20'))
plt.add_trace(go.Scatter(x=x, y=[prob(i, N = 50) for i in x], name = 'N = 50'))
plt.update_layout(title='Pr(helix/random) over Temperature by residue count',
                  yaxis = dict(range=[0,2]), yaxis_title = 'Probability(helix/random)',
                  xaxis = dict(range=[320-100, 320+100]), xaxis_title = 'Temperature')
plt.show()

# Changing Enthalpy
plt_ep = go.Figure()
plt_ep.add_trace(go.Scatter(x=x, y=[prob(i, ep = 1.15e-21) for i in x], name='Lower Enthalpy'))
plt_ep.add_trace(go.Scatter(x=x, y=[prob(i) for i in x], name='Control'))
plt_ep.add_trace(go.Scatter(x=x, y=[prob(i, ep = 2e-21) for i in x], name='Higher Enthalpy'))
plt_ep.update_layout(title='Pr(helix/random) over Temperature by changing enthalpy of H Bonds',
                     yaxis = dict(range=[0, 2]), yaxis_title='Probability(helix/random)',
                     xaxis = dict(range=[320-100, 320+100]), xaxis_title = 'Temperature')
plt_ep.show()

# Changing Entropy
plt_s = go.Figure()
plt_s.add_trace(go.Scatter(x=x, y=[prob(i, s = 4e-24) for i in x], name = 'Lower Entropy'))
plt_s.add_trace(go.Scatter(x=x, y=[prob(i) for i in x], name = 'Control'))
plt_s.add_trace(go.Scatter(x=x, y=[prob(i, s = 6.1e-24) for i in x], name = 'Higher Entropy'))
plt_s.update_layout(title='Pr(helix/random) over Temperature by changing entropy',
                    yaxis = dict(range=[0,2]), yaxis_title = 'Probability(helix/random)',
                    xaxis = dict(range=[320-100, 320+100]), xaxis_title = 'Temperature')
plt_s.show()
