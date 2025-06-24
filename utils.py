'''
Utilities
'''

# Imports
import sciris as sc
import hpvsim as hpv


def set_font(size=None, font='Libertinus Sans'):
    ''' Set a custom font '''
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def shrink_calib(calib, n_results=100):
    plot_indices = calib.df.iloc[0:n_results, 0].values
    calib.analyzer_results = [calib.analyzer_results[i] for i in plot_indices]
    calib.sim_results = [calib.sim_results[i] for i in plot_indices]
    calib.extra_sim_results = [calib.extra_sim_results[i] for i in plot_indices]
    calib.target_data = calib.target_data
    calib.df = calib.df.iloc[0:n_results, ]
    return calib

 