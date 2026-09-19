"""Table S3: threshold TxV per-dose price by scenario, both comparators.

Reads results/scens_long.csv.gz (per scenario x sim x year x metric) and
solves for the TxV per-dose price P* at which each TxV-carrying scenario
matches a comparator on cost-effectiveness. Two comparators (status-quo
S&T&T 18% and best non-TxV S&T 70%), three WTP tiers per DALY.

Costs are applied to intervention counts (screens, ablations, LEEPs,
radiation courses, prophylactic doses); TxV dose price P is the unknown.
All streams are discounted at 3%/yr from 2030.
"""
import numpy as np
import pandas as pd


CUM_START_YEAR = 2030
DISCOUNT = 0.03

TXV_SCENARIOS = [
    'S&TxV 18%',           'S&TxV 35%',           'S&TxV 70%',
    'S&TxV&T&T 18%',       'S&TxV&T&T 35%',       'S&TxV&T&T 70%',
    'Mass TxV 90/0, 18%',  'Mass TxV 90/0, 35%',  'Mass TxV 90/0, 70%',
    'Mass TxV 50/90, 18%', 'Mass TxV 50/90, 35%', 'Mass TxV 50/90, 70%',
]

COMPARATORS = ['S&T&T 18%', 'S&T 70%']

# Unit costs (USD 2024). base = point estimate; low/high = literature range.
# Sources: Campos et al 2020 (Kenya/Uganda CEA), Sy et al 2022 (Rwanda-adjacent),
# WHO-CHOICE East Africa service costs, GAVI HPV vaccine pricing.
UNIT_COST = {
    'screens':            dict(base=10,   low=5,   high=15),
    'ablations':          dict(base=25,   low=15,  high=40),
    'leeps':              dict(base=100,  low=50,  high=150),
    'cancer_treatments':  dict(base=1500, low=500, high=3000),
    'vaccinations':       dict(base=7,    low=5,   high=10),
}

# WTP per DALY averted (USD). Opportunity-cost anchor: Ochalek et al 2018
# sub-Saharan Africa estimates. GDP anchors: Rwanda GDP per capita ~$900
# (World Bank 2024).
WTP_DALY = [('opp_cost_130', 130), ('gdp_half_450', 450), ('gdp_full_900', 900)]


def _discount_factor(year, base=CUM_START_YEAR, rate=DISCOUNT):
    return 1.0 / (1.0 + rate) ** (np.asarray(year) - base)


def _discounted_totals(long):
    df = long[long.year >= CUM_START_YEAR].copy()
    df['dv'] = df['value'] * _discount_factor(df['year'].values)
    return (df.groupby(['scenario', 'sim', 'metric'])['dv'].sum()
              .rename('disc_total').reset_index())


def _wide(totals):
    return (totals.set_index(['scenario', 'sim', 'metric'])['disc_total']
                  .unstack('metric').fillna(0.0))


def _cost_excl_txv(wide, field='base'):
    cost = np.zeros(len(wide))
    for metric, prices in UNIT_COST.items():
        if metric in wide.columns:
            cost += wide[metric].values * prices[field]
    return pd.Series(cost, index=wide.index)


def _paired_stats(scen_series, comp_series, sign):
    """Median + 10-90 quantiles of paired diff sign*(a - b) per sim."""
    diff = sign * (scen_series - comp_series)
    return dict(med=diff.median(), lo=diff.quantile(0.10), hi=diff.quantile(0.90))


def _p_star(d_cost_excl, d_effect, d_doses, wtp):
    if d_doses <= 0:
        return np.nan
    return (wtp * d_effect - d_cost_excl) / d_doses


def build(resfolder='results'):
    long = pd.read_csv(f'{resfolder}/scens_long.csv.gz')
    wide = _wide(_discounted_totals(long))

    for m in list(UNIT_COST.keys()) + ['txvs', 'cancers', 'dalys']:
        if m not in wide.columns:
            wide[m] = 0.0

    # Precompute cost_excl per (scenario, sim) under each unit-cost field.
    cost_excl = {f: _cost_excl_txv(wide, f) for f in ('base', 'low', 'high')}

    rows = []
    for scen in TXV_SCENARIOS:
        s_lev = wide.index.get_level_values('scenario') == scen
        s_sub = wide.loc[s_lev].reset_index('scenario', drop=True)
        for comp in COMPARATORS:
            c_lev = wide.index.get_level_values('scenario') == comp
            c_sub = wide.loc[c_lev].reset_index('scenario', drop=True)

            d_txv_stats  = _paired_stats(s_sub['txvs'],    c_sub['txvs'],    +1)
            d_daly_stats = _paired_stats(c_sub['dalys'],   s_sub['dalys'],   +1)  # averted
            d_canc_stats = _paired_stats(c_sub['cancers'], s_sub['cancers'], +1)  # averted

            row = dict(
                scenario=scen,
                comparator=comp,
                dalys_averted_med=d_daly_stats['med'],
                dalys_averted_lo=d_daly_stats['lo'],
                dalys_averted_hi=d_daly_stats['hi'],
                cancers_averted_med=d_canc_stats['med'],
                cancers_averted_lo=d_canc_stats['lo'],
                cancers_averted_hi=d_canc_stats['hi'],
                txv_doses_med=d_txv_stats['med'],
            )
            # ICER at P=0: cost per DALY averted with free TxV.
            base_cost_diff = (cost_excl['base'].loc[s_lev].values.mean()
                              - cost_excl['base'].loc[c_lev].values.mean())
            if d_daly_stats['med'] > 0:
                row['icer_p0'] = base_cost_diff / d_daly_stats['med']
            else:
                row['icer_p0'] = np.nan

            for label, wtp in WTP_DALY:
                # Base + range across cost inputs.
                p_vals = {}
                for field in ('base', 'low', 'high'):
                    cs = cost_excl[field]
                    diff_cost = (cs.loc[s_lev].values.mean()
                                 - cs.loc[c_lev].values.mean())
                    p_vals[field] = _p_star(diff_cost, d_daly_stats['med'],
                                            d_txv_stats['med'], wtp)
                row[f'P_{label}_base'] = p_vals['base']
                row[f'P_{label}_low']  = min(p_vals.values())
                row[f'P_{label}_high'] = max(p_vals.values())

            rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(f'{resfolder}/table_s3_costing.csv', index=False, float_format='%.2f')
    print(f'Wrote {resfolder}/table_s3_costing.csv ({len(out)} rows)')
    return out


if __name__ == '__main__':
    build()
