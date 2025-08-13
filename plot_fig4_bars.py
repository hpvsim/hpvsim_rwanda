"""
Plot residual burden of cervical cancer in Rwanda under vaccination scenarios
"""
import pandas as pd


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig4(end_year=2060):
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    ut.set_font(20)
    fig = pl.figure(figsize=(18, 9))
    gs1 = pl.GridSpec(1, 3, left=0.05, right=0.99, bottom=0.60, top=0.95, wspace=0)
    gs2 = pl.GridSpec(1, 3, left=0.05, right=0.99, bottom=0.12, top=0.45, wspace=0)

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    fi = sc.findinds(mbase.year, 2025)[0]
    ei = sc.findinds(mbase.year, end_year)[0]+1

    # Remove some results
    msim_dict.pop('No interventions')
    msim_dict.pop('Baseline')
    delkeys = [k for k in msim_dict.keys() if 'Male' in k or 'HIV' in k]
    for k in delkeys:
        msim_dict.pop(k)

    # Assemble cumulative results
    cum_res = {k:dict() for k in ['18', '35', '70']}
    label_map = {
        'Baseline': 'Status quo',
        'TxV 50/90 with 18% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 18% screening': 'Virus-\nclearing\ntherapeutic',
        'TxV 50/90 with 35% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 35% screening': 'Virus-\nclearing\ntherapeutic',
        'TxV 50/90 with 70% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 70% screening': 'Virus-\nclearing\ntherapeutic',
        'Mass vx 18%': 'One-time\nvax',
        'S&T 18%': 'Screen\n& treat',
        'Mass vx 35%': 'One-time\nvax',
        'S&T 35%': 'Screen\n& treat',
        'Mass vx 70%': 'One-time\nvax',
        'S&T 70%': 'Screen\n& treat',
    }
    resname = 'cancers'
    for sname, scen in msim_dict.items():
        sl = sname.split('%')[0][-2:]
        cum_res[sl][sname] = scen[resname].values[fi:ei].sum()

    new_order = [
        'One-time\nvax',
        'Virus-\nclearing\ntherapeutic',
        'Lesion-\nregressing\ntherapeutic',
        'Screen\n& treat'
     ]
    # Make a bar plot
    perc_averted = dict()
    for pn, slevel in enumerate(['18', '35', '70']):
        ax = fig.add_subplot(gs1[pn])

        # Rearrage order of dict items
        remap = {label_map[k]:v for k,v in cum_res[slevel].items()}
        remap = {k:remap[k] for k in new_order if k in remap}
        bars = list(remap.values())
        x = np.arange(len(remap))
        ax.bar(x, bars, color=sc.gridcolors(len(bars)))
        ax.set_xticks(x)
        ax.set_xticklabels(remap.keys(), ha='center')
        ax.set_title('Cumulative cancers')
        sc.SIticks(ax)
        ax.set_ylim(bottom=0, top=80_000)
        if pn != 0: ax.set_yticklabels([])
        text_loc = 80_000*.95
        ax.text(1.5, text_loc, f'{slevel}% coverage', va="center", ha="center")

        ax = fig.add_subplot(gs2[pn])
        new_bars = [mbase['cancers'].values[fi:ei].sum()-b for b in bars]
        perc_averted[slevel] = [100*b/mbase['cancers'].values[fi:ei].sum() for b in new_bars]
        ax.bar(x, new_bars, color=sc.gridcolors(len(new_bars)))
        ax.set_xticks(x)
        ax.set_xticklabels(remap.keys(), ha='center')
        ax.set_title('Cancers averted')
        sc.SIticks(ax)
        ax.set_ylim(bottom=0, top=30_000)
        if pn != 0: ax.set_yticklabels([])
        text_loc = 30_000*.95
        ax.text(1.5, text_loc, f'{slevel}% coverage', va="center", ha="center")

    fig.tight_layout()
    fig_name = f'figures/fig4_bars_{end_year}.png'
    sc.savefig(fig_name, dpi=100)

    return perc_averted, new_bars


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    perc_averted, num_averted = plot_fig4(end_year=2060)




