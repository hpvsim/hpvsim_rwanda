import numpy as np, pandas as pd
def load(b):
    return (pd.read_csv(f'results/{b}/scens_timeseries.csv'),
            pd.read_csv(f'results/{b}/scens_cumulative.csv'))
v2ts, v2c = load('v2.3.1_baseline'); v3ts, v3c = load('v3.0_baseline')
def val(ts, scen, metric, yr):
    s = ts[(ts.scenario == scen) & (ts.metric == metric)]
    if len(s) == 0: return np.nan
    i = (s.year - yr).abs().idxmin(); return float(s.loc[i, 'value'])
scens = ['No interventions', 'Baseline', 'S&T&T 70%', 'S&T 70%']
print('=== ASR cancer incidence (per 100k, WHO-std) ===')
print(f'{"scenario":18} {"yr":>4} {"v2.3.1":>8} {"v3.0":>8}')
for sc in scens:
    for yr in [2020, 2035, 2050]:
        print(f'{sc:18} {yr:>4} {val(v2ts,sc,"asr_cancer_incidence",yr):>8.1f} {val(v3ts,sc,"asr_cancer_incidence",yr):>8.1f}')
print('\n=== HIV+ vs HIV- crude cancer incidence (per 100k), Baseline ===')
for yr in [2015, 2020, 2035, 2050]:
    print(f'{yr}: v2 HIV+={val(v2ts,"Baseline","cancer_incidence_with_hiv",yr):.0f} '
          f'HIV-={val(v2ts,"Baseline","cancer_incidence_no_hiv",yr):.0f} | '
          f'v3 HIV+={val(v3ts,"Baseline","cancer_incidence_with_hiv",yr):.0f} '
          f'HIV-={val(v3ts,"Baseline","cancer_incidence_no_hiv",yr):.0f}')
print('\n=== HIV rate ratio (HIV+/HIV-) mean 2015-2035 ===')
for eng, ts in [('v2.3.1', v2ts), ('v3.0', v3ts)]:
    hp = ts[(ts.scenario=='Baseline')&(ts.metric=='cancer_incidence_with_hiv')&(ts.year>=2015)&(ts.year<=2035)].value.mean()
    hn = ts[(ts.scenario=='Baseline')&(ts.metric=='cancer_incidence_no_hiv')&(ts.year>=2015)&(ts.year<=2035)].value.mean()
    print(f'{eng}: HIV+ mean={hp:.0f} HIV- mean={hn:.0f} RR={hp/hn if hn else float("nan"):.1f}')
print('\n=== cumulative cancers 2025-2050 relative to No-interventions (scale-invariant) ===')
def cum(c, sc):
    r = c[(c.scenario == sc) & (c.metric == 'cancers')]; return float(r.value.iloc[0])
for eng, c in [('v2.3.1', v2c), ('v3.0', v3c)]:
    base = cum(c, 'No interventions')
    print(eng, {sc: round(cum(c, sc) / base, 3) for sc in scens})
print('\n=== mean ASR 2016-2050 per scenario ===')
for sc in scens:
    m2 = v2ts[(v2ts.scenario==sc)&(v2ts.metric=='asr_cancer_incidence')&(v2ts.year>=2016)&(v2ts.year<=2050)].value.mean()
    m3 = v3ts[(v3ts.scenario==sc)&(v3ts.metric=='asr_cancer_incidence')&(v3ts.year>=2016)&(v3ts.year<=2050)].value.mean()
    print(f'{sc:18} v2={m2:6.1f} v3={m3:6.1f}')
