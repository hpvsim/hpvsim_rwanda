"""Build Table S2: cumulative intervention counts per scenario 2030-2100 (R2.7).

Reads results/scens_cumulative.csv (produced by run_scenarios.py) and writes
results/table_s2.csv.

Columns: scenario, cumulative cancers, screens, ablations, LEEP procedures,
radiation courses, therapeutic vaccine doses, prophylactic vaccine doses.

Campaign-era 'excisions' from make_st_older are excluded pending resolution of
the older-cohort counter accumulation issue (rolling to-do).
"""
import pandas as pd

from utils import load_scens, get_cum


SCEN_ORDER = [
    'No interventions',
    'S&T&T 18%',
    'S&T&T 35%',
    'S&T&T 70%',
    'S&T 18%',
    'S&T 35%',
    'S&T 70%',
    'S&T 70%, 50% LTFU',
    'S&TxV&T&T 18%',
    'S&TxV&T&T 35%',
    'S&TxV&T&T 70%',
    'S&TxV 18%',
    'S&TxV 35%',
    'S&TxV 70%',
    'HPV-Faster 18%',
    'HPV-Faster 35%',
    'HPV-Faster 70%',
    'Mass TxV 90/0, 18%',
    'Mass TxV 90/0, 35%',
    'Mass TxV 90/0, 70%',
    'Mass TxV 50/90, 18%',
    'Mass TxV 50/90, 35%',
    'Mass TxV 50/90, 70%',
]

METRICS = [
    ('cancers',            'cumulative_cancers'),
    ('screens',            'screens'),
    ('ablations',          'ablations'),
    ('leeps',              'LEEP_procedures'),
    ('cancer_treatments',  'radiation_courses'),
    ('txvs',               'therapeutic_vaccine_doses'),
    ('vaccinations',       'prophylactic_vaccine_doses'),
]


def main():
    _, cum, _, _ = load_scens('results')
    rows = []
    for s in SCEN_ORDER:
        row = {'scenario': s}
        for src, dst in METRICS:
            try:
                med, _, _ = get_cum(cum, s, src)
                row[dst] = med
            except Exception:
                row[dst] = 0.0
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv('results/table_s2.csv', index=False, float_format='%.0f')
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
