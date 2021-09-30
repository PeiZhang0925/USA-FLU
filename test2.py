import pandas as pd
import numpy as np

import matplotlib
from matplotlib import pyplot as plt

from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr

data_all = pd.read_csv("/Users/pei/Dropbox (linweit)/ehrg/users/Pei/PJT/USA_FLU/codes/usaflu_raw.csv")


# Function for preparing the values of the variables specified by 'var_names'
# in the state specified by 'state' as numpy array of shape (T, N), where
# T is the number of time steps and N the number of variables
def get_data(data_all, state, var_names=['fluP', 'o3', 'season', 'temp', 'ah']):
    # Select the state
    data_out = data_all.loc[data_all["state"] == state]

    # Select the columns
    data_out = data_out[var_names]

    # Turn into numpy array
    data_out = data_out.values

    # Return
    return data_out


# Function for preparing the dictionary that will be passed to the 'selected_links'
# argument of run_pcmciplus() and run_pcmci():
# 'fluP' can be caused by all variables
# 'o3' can be caused by all variables but 'fluP'
# other variables can only be caused by themselves
def get_selected_links(var_names):
    # Build dictionary
    selected_links = {}
    for idx, var in enumerate(var_names):
        if var == 'fluP':
            selected_links[idx] = [(other_idx, -tau) for other_idx, _ in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1)]
        elif var == 'o3':
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) if other_var != 'fluP']
        else:
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

    # Return
    return selected_links


# Specify the state and analyzed variables
state = "New York"
var_names = ['fluP', 'o3', 'season', 'temp', 'ah']

# Prepare the tigramite dataframe
data = get_data(data_all, state=state, var_names=var_names)
dataframe = pp.DataFrame(data, var_names=var_names, missing_flag=999.)

# Instantiate a PCMCI oject using the ParCorr conditional independence test
parcorr = ParCorr()
pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, verbosity=0)

# Specify parameters of run_pcmciplus
tau_min = 0
tau_max = 2
pc_alpha = None  # pc_alpha = [0.01, 0.1, 0.2, 0.3], pc_alpha = 0.1
selected_links = get_selected_links(var_names)

# Run PCMCI^+ with these parameters
results = pcmci.run_pcmciplus(tau_min=tau_min,
                              tau_max=tau_max,
                              pc_alpha=pc_alpha,
                              selected_links=selected_links)

# Plot the summary graph (not resolved in time)
# NOTE:
# When working with PCMCI^+, pass the estimated graph (here stored in
# results['graph']) directly to the 'link_matrix' argument of the
# plotting functions. Do not use return_significant_links.
tp.plot_graph(link_matrix=results['graph'],
              val_matrix=results['val_matrix'],
              var_names=var_names)
plt.show()

# Plot the time series graph (resolved in time)
tp.plot_time_series_graph(link_matrix=results['graph'],
                          val_matrix=results['val_matrix'],
                          var_names=var_names)
plt.show()