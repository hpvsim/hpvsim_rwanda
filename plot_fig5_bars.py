"""
Plot residual burden of cervical cancer in Rwanda under vaccination scenarios
"""
import pandas as pd


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig4(end_year=2100):
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    ut.set_font(20)
    fig, axes = pl.subplots(nrows=2, ncols=1, figsize=(12, 10))
    axes = axes.ravel()
    # gs1 = pl.GridSpec(2, 3, left=0.05, right=0.99, bottom=0.60, top=0.95, wspace=0)
    # gs2 = pl.GridSpec(2, 3, left=0.05, right=0.99, bottom=0.12, top=0.45, wspace=0)

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    fi = sc.findinds(mbase.year, 2025)[0]
    ei = sc.findinds(mbase.year, end_year)[0]+1
    ei2 = sc.findinds(mbase.year, 2060)[0]+1

    # Remove some results
    msim_dict.pop('No interventions')
    msim_dict.pop('Baseline')
    delkeys = [k for k in msim_dict.keys() if 'Male' in k or 'HIV' in k]
    for k in delkeys:
        msim_dict.pop(k)

    # Assemble cumulative results
    cum_res = {k:dict() for k in ['18', '35', '70']}
    cum_res_lb = {k:dict() for k in ['18', '35', '70']}
    cum_res_ub = {k:dict() for k in ['18', '35', '70']}

    # Map labels
    vclab = 'Virus-clearing'
    lrlab = 'Lesion-regressing'

    label_map = {
        'Baseline': 'Status quo',
        'TxV 50/90 within screening, 18%': 'Screen & TxV\n'+lrlab,
        'TxV 90/50 within screening, 18%': 'Screen & TxV\n'+vclab,
        'TxV 50/90 within screening, 35%': 'Screen & TxV\n'+lrlab,
        'TxV 90/50 within screening, 35%': 'Screen & TxV\n'+vclab,
        'TxV 50/90 within screening, 70%': 'Screen & TxV\n'+lrlab,
        'TxV 90/50 within screening, 70%': 'Screen & TxV\n'+vclab,
        'Mass TxV 90/50, 18%': 'Mass TxV\n'+vclab,
        'Mass TxV 50/90, 18%': 'Mass TxV\n'+lrlab,
        'Mass TxV 90/50, 35%': 'Mass TxV\n'+vclab,
        'Mass TxV 50/90, 35%': 'Mass TxV\n'+lrlab,
        'Mass TxV 90/50, 70%': 'Mass TxV\n'+vclab,
        'Mass TxV 50/90, 70%': 'Mass TxV\n'+lrlab,
        'HPV-Faster 18%': 'HPV-Faster',
        'HPV-Faster 35%': 'HPV-Faster',
        'HPV-Faster 70%': 'HPV-Faster',
        # 'S&T 18%': 'Screen\n& treat',
        # 'S&T 35%': 'Screen\n& treat',
        # 'S&T 70%': 'Screen\n& treat',
    }
    resname = 'cancers'
    for sname, scen in msim_dict.items():
        sl = sname.split('%')[0][-2:]
        cum_res[sl][sname] = scen[resname].values[fi:ei].sum()
        cum_res_lb[sl][sname] = scen[resname].low[fi:ei].sum()
        cum_res_ub[sl][sname] = scen[resname].high[fi:ei].sum()

    new_order = [
        # 'Status quo',
        'HPV-Faster',
        'Mass TxV\n'+vclab,
        'Mass TxV\n'+lrlab,
        'Screen & TxV\n'+ vclab,
        'Screen & TxV\n'+ lrlab,
     ]

    x = np.arange(len(new_order))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0
    multiplier2 = 0
    colors = sc.gridcolors(3)

    # Make bar plots
    perc_averted = dict()
    for pn, slevel in enumerate(['18', '35', '70']):
        ax = axes[0]

        offset = width * multiplier

        # Rearrage order of dict items
        remap = {label_map[k]:v for k,v in cum_res[slevel].items()}
        remap = {k:remap[k] for k in new_order if k in remap}
        remap_lb = {label_map[k]:v for k,v in cum_res_lb[slevel].items()}
        remap_lb = {k:remap_lb[k] for k in new_order if k in remap}
        remap_ub = {label_map[k]:v for k,v in cum_res_ub[slevel].items()}
        remap_ub = {k:remap_ub[k] for k in new_order if k in remap}
        # Calculate yerr
        yerr = np.array([remap_ub[k]-remap_lb[k] for k in remap.keys()])

        bars = list(remap.values())
        ax.bar(x + offset-width, bars, width, color=colors[pn])  #, yerr=yerr)
        # Add y error bars

        ax.set_xticks(x)
        ax.set_xticklabels(remap.keys(), ha='center')
        ax.set_title('Cumulative cancers 2025-2100')
        sc.SIticks(ax)
        top = 90_000
        ax.set_ylim(bottom=0, top=top)

        # Add gridlines
        ax.grid(axis='y', linestyle='--', alpha=0.5)
        # Add bar labels
        for i, v in enumerate(bars):
            ax.text(i + offset - width, v + top * 0.01, f'{v:,.0f}', ha='center', va='bottom', fontsize=14)

        ax = axes[1]
        offset2 = width * multiplier2
        new_bars = [mbase['cancers'].values[fi:ei].sum()-b for b in bars]
        perc_averted[slevel] = [100*b/mbase['cancers'].values[fi:ei].sum() for b in new_bars]
        ax.bar(x+offset2-width, new_bars, width, color=colors[pn], label=f'{slevel}% coverage')
        ax.set_xticks(x)
        ax.set_xticklabels(remap.keys(), ha='center')
        ax.set_title('Cancers averted 2025-2100')
        sc.SIticks(ax)
        top = 50_000
        ax.set_ylim(bottom=0, top=top)
        ax.legend(loc='upper left', frameon=False)

        # Add gridlines
        ax.grid(axis='y', linestyle='--', alpha=0.5)
        # Add bar labels
        for i, v in enumerate(new_bars):
            ax.text(i + offset2 - width, v + top * 0.01, f'{v:,.0f}', ha='center', va='bottom', fontsize=14)

        multiplier += 1
        multiplier2 += 1

    fig.tight_layout()
    fig_name = f'figures/fig5_bars.png'
    sc.savefig(fig_name, dpi=100)

    return perc_averted, new_bars


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    perc_averted, num_averted = plot_fig4()




