# Imports

import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()
import pyarrow.feather as feather

import matplotlib
from matplotlib import pyplot as plt
## %matplotlib inline
## use `%matplotlib notebook` for interactive figures
# plt.style.use('ggplot')
import sklearn

import tigramite
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb


def my_parcorr(pc_sig,final_sig,tau_min,tau_max,var_names,i):
    var_names=var_names
    dat = dat0
    dat1=dat[var_names]
    dat1.rename(columns={'fluP': r'$Flu$',
                          'o3': r'$O_{3}$',
                          'season': r'$Season$',
                         'temp': r'$Temp.$',
                         'ah': r'$Humid.$'}, inplace=True)
    var_names=[r'$Flu$',r'$O_{3}$',r'$Season$',
    r'$Temp.$',r'$Humid.$']
    dat2 = np.array(dat1)

    dataframe = pp.DataFrame(dat2,
                         datatime = np.arange(len(dat2)) ,
                         missing_flag=999.,
                         var_names=var_names)
    parcorr = ParCorr(significance='analytic')
    pcmci = PCMCI(
        dataframe=dataframe,
        cond_ind_test=parcorr,
        verbosity=0)
    tau_max = tau_max
    tau_min = tau_min
    pc_alpha = pc_sig
    alpha_level = final_sig
    pcmci.verbosity = 0

    selected_links = {}
    for j in range(len(var_names)):
        if j in [var_names.index(r'$Season$')]:
            selected_links[j] = [(var, -lag) for var in [var_names.index(r'$Season$')]
                         for lag in range(tau_min, tau_max + 1)]
        elif j in [var_names.index(r'$Temp.$')]:
            selected_links[j] = [(var, -lag) for var in [var_names.index(r'$Temp.$')]
                         for lag in range(tau_min, tau_max + 1)]
        elif j in [var_names.index(r'$Humid.$')]:
            selected_links[j] = [(var, -lag) for var in [var_names.index(r'$Humid.$')]
                         for lag in range(tau_min, tau_max + 1)]
        # elif j in [var_names.index("o3")]:
        #     selected_links[j] = [(var, -lag) for var in [var_names.index("o3")]
        #                  for lag in range(tau_min, tau_max + 1)]
        elif j in [0]:
            selected_links[j] = [(var, -lag) for var in range(len(var_names))
                         for lag in range(tau_min, tau_max + 1)]
        else:
            selected_links[j] = [(var, -lag) for var in [i for i in range(len(var_names)) if i != 0 ]
                         for lag in range(tau_min, tau_max + 1)]

    results = pcmci.run_pcmciplus(tau_min=0, tau_max=tau_max, pc_alpha=pc_alpha, selected_links=selected_links)

    link_matrix = pcmci.return_significant_links(pq_matrix=results['p_matrix'],
                        val_matrix=results['val_matrix'], alpha_level=alpha_level)['link_matrix']

    c3d=results['val_matrix']
    c2d=c3d.transpose(2,0,1).reshape(tau_max+1,-1)
    c_df=pd.DataFrame(c2d)
    c_df_t=c_df.transpose()
    c_df_t["state"]=i
    feather.write_feather(c_df_t, 'c_df%s.feather' %i)

    p3d=results['p_matrix']
    p2d=p3d.transpose(2,0,1).reshape(tau_max+1,-1)
    p_df=pd.DataFrame(p2d)
    p_df_t=p_df.transpose()
    p_df_t["state"]=i
    feather.write_feather(p_df_t, 'p_df%s.feather' %i)

    l3d=link_matrix
    l2d=l3d.transpose(2,0,1).reshape(tau_max+1,-1)
    l_df=pd.DataFrame(l2d)
    l_df_t=l_df.transpose()
    l_df_t["state"]=i
    feather.write_feather(l_df_t, 'l_df%s.feather' %i)

    tp.plot_graph(
        arrow_linewidth=12.0,
        figsize = (10,5),
        vmin_edges=-0.5,
        vmax_edges=0.5,
        node_label_size = 15,
        link_label_fontsize = 10,
        val_matrix=results['val_matrix'],
        link_matrix=link_matrix,
        var_names=var_names,
        link_colorbar_label='cross-MCI (edges)',
        node_colorbar_label='auto-MCI (nodes)',
        label_fontsize=15,
        network_lower_bound=0.2,
        show_colorbar=1
        );
    plt.suptitle("Overall",size=18,weight="semibold",verticalalignment="top",horizontalalignment="center")
    plt.savefig('/Users/pei/Dropbox (linweit)/ehrg/users/Pei/PJT/USA_FLU/codes/PCMCIres_nation.png',dpi=500)
    plt.show()

dat0 = pd.read_csv("/Users/pei/Dropbox (linweit)/ehrg/users/Pei/PJT/USA_FLU/codes/usaflu_norm.csv", header=0, index_col=0)

states = ["OVERALL"]

for i in states:
  my_parcorr(None,0.05,0,2,var_names=['fluP','o3','season','temp','ah'],i=i)