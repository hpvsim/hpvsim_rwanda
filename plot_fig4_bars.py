"""
Plot residual burden of cervical cancer in Rwanda under vaccination scenarios
"""
import pandas as pd


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig4():
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    ut.set_font(20)
    fig = pl.figure(figsize=(15, 9))
    gs1 = pl.GridSpec(1, 3, left=0.05, right=0.99, bottom=0.60, top=0.95, wspace=0)
    gs2 = pl.GridSpec(1, 3, left=0.05, right=0.99, bottom=0.12, top=0.45, wspace=0)

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    fi = sc.findinds(mbase.year, 2025)[0]

    # Remove some results
    msim_dict.pop('No interventions')
    msim_dict.pop('Baseline')
    delkeys = [k for k in msim_dict.keys() if 'Male' in k or 'HIV' in k]
    for k in delkeys:
        msim_dict.pop(k)

    # Assemble cumulative results
    cum_res = {k:dict() for k in ['10', '50', '90']}
    label_map = {
        'Baseline': 'Status quo',
        'TxV 50/90 with 10% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 10% screening': 'Virus-\nclearing\ntherapeutic',
        'TxV 50/90 with 50% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 50% screening': 'Virus-\nclearing\ntherapeutic',
        'TxV 50/90 with 90% screening': 'Lesion-\nregressing\ntherapeutic',
        'TxV 90/50 with 90% screening': 'Virus-\nclearing\ntherapeutic',
        'Mass vx 10%': 'One-time\nvax',
        'S&T 10%': 'Screen\n& treat',
        'Mass vx 50%': 'One-time\nvax',
        'S&T 50%': 'Screen\n& treat',
        'Mass vx 90%': 'One-time\nvax',
        'S&T 90%': 'Screen\n& treat',
    }
    resname = 'cancers'
    for sname, scen in msim_dict.items():
        sl = sname.split('%')[0][-2:]
        cum_res[sl][sname] = scen[resname].values[fi:].sum()

    # Make a bar plot
    for pn, slevel in enumerate(['10', '50', '90']):
        ax = fig.add_subplot(gs1[pn])

        # Rearrage order of dict items
        sorted_res = sc.odict(sorted(cum_res[slevel].items(), key=lambda x: x[1], reverse=True))
        bars = list(sorted_res.values())
        labels = [label_map[k] for k in sorted_res.keys()]
        x = np.arange(len(labels))
        ax.bar(x, bars, color=sc.gridcolors(len(bars)))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, ha='center')
        ax.set_title('Cumulative cancers')
        sc.SIticks(ax)
        ax.set_ylim(bottom=0, top=80_000)
        if pn != 0: ax.set_yticklabels([])
        text_loc = 80_000*.95
        ax.text(1.5, text_loc, f'{slevel}% coverage', va="center", ha="center")

        ax = fig.add_subplot(gs2[pn])
        new_bars = [mbase['cancers'].values[fi:].sum()-b for b in bars]
        perc_averted = [100*b/mbase['cancers'].values[fi:].sum() for b in new_bars]
        ax.bar(x, new_bars, color=sc.gridcolors(len(new_bars)))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, ha='center')
        ax.set_title('Cancers averted')
        sc.SIticks(ax)
        ax.set_ylim(bottom=0, top=50_000)
        if pn != 0: ax.set_yticklabels([])
        text_loc = 50_000*.95
        ax.text(1.5, text_loc, f'{slevel}% coverage', va="center", ha="center")

    fig.tight_layout()
    fig_name = 'figures/fig4_bars.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig4()




