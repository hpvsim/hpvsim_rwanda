"""Test hypothesis: v2 used screen_change_year=2025, v3 uses 2027. Does moving
back to 2025 restore the 9K averted-cancer delta between S&T&T 18% and 70%?

Runs S&T&T 18%, S&T&T 70%, S&T&T 70% (2025 change year) on 3 reps and prints
the cumulative cancer window sums.
"""
import numpy as np
import pandas as pd
import sciris as sc
import run_scenarios as rsc
import run_sim as rs
import interventions as intv


def run_one(name, top_par, change_year=2027, fut_cov=0.18):
    intvs = intv.make_st(future_screen_cov=fut_cov,
                          screen_change_year=change_year,
                          end_year=2100)
    sim = rs.make_sim(add_st=False, interventions=intvs,
                      stop=2100, calib_pars=top_par)
    sim.run(verbose=0)
    r = sim.results.all_hpv.new_cancers
    tv = np.asarray(r.timevec.years if hasattr(r.timevec, 'years') else r.timevec)
    mask = (tv >= 2025) & (tv < 2100)
    return float(np.asarray(r.values)[mask].sum())


if __name__ == '__main__':
    top_pars = rsc._top_pars(3)
    print('\n=== Q1 test: does screen_change_year shift the S&T&T scale-up delta? ===')
    rows = []
    for rep, tp in enumerate(top_pars):
        print(f'\nrep {rep}')
        T = sc.timer()
        # S&T&T 18% (change year irrelevant since same coverage)
        c18 = run_one('sTT18', dict(tp), fut_cov=0.18)
        T.toc('sTT18')
        # S&T&T 70% at 2027 change year (current)
        c70_2027 = run_one('sTT70_2027', dict(tp), change_year=2027, fut_cov=0.70)
        T.toc('sTT70_2027')
        # S&T&T 70% at 2025 change year (MS's v2)
        c70_2025 = run_one('sTT70_2025', dict(tp), change_year=2025, fut_cov=0.70)
        T.toc('sTT70_2025')
        # S&T&T 70% at 2020 change year (immediate 70%)
        c70_2020 = run_one('sTT70_2020', dict(tp), change_year=2020, fut_cov=0.70)
        T.toc('sTT70_2020')
        rows.append(dict(rep=rep, c18=c18, c70_2027=c70_2027,
                          c70_2025=c70_2025, c70_2020=c70_2020,
                          avert_2027=c18-c70_2027, avert_2025=c18-c70_2025,
                          avert_2020=c18-c70_2020))
        print(f'  18%:                        {c18:,.0f}')
        print(f'  70% change_year=2027:       {c70_2027:,.0f}  (avert {c18-c70_2027:,.0f})')
        print(f'  70% change_year=2025:       {c70_2025:,.0f}  (avert {c18-c70_2025:,.0f})')
        print(f'  70% change_year=2020:       {c70_2020:,.0f}  (avert {c18-c70_2020:,.0f})')
    df = pd.DataFrame(rows)
    df.to_csv('results/diagnostic/q1_change_year.csv', index=False)
    print('\n--- Medians ---')
    print(df.median().to_string())
