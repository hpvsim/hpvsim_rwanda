"""Light re-fit of ONLY the HIV-differential params (rel_sev_gt200, and
optionally rel_sev_lt200) in RWANDA_HIV_EFFECTS, at ms_agent_ratio=5.

NOT a full Optuna: a 1-D grid over rel_sev_gt200 (the dominant lever per project
history) holding the HPV natural-history scalars and rel_sus fixed. Uses the
same pooled HIV+ estimator and 2010-2019 window as the canonical calibration
(sum weighted cancers / sum weighted female-years), scale-weighted for grow
multiscale. Targets: HIV+ agg 33.0, HIV- agg 13.1.

Usage:
    .venv/Scripts/python.exe refit_hiv_effects.py [n_agents] [n_seeds] [grid_csv]
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')
import sys
import copy

import numpy as np
import sciris as sc

import run_sim as rs
import verify_hiv_fullscale as vf

N_AGENTS = int(sys.argv[1]) if len(sys.argv) > 1 else 15_000
N_SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
GRID = ([float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3
        else [1.4, 1.6, 1.8, 1.98, 2.2])   # rel_sev_gt200 candidates

TARGET_HIVP = 33.0
TARGET_HIVN = 13.1


def _effects_with(gt200):
    e = copy.deepcopy(rs.rc.RWANDA_HIV_EFFECTS)
    e['rel_sev'] = dict(e['rel_sev'])
    e['rel_sev']['gt200'] = gt200
    return e


if __name__ == '__main__':
    import hpvsim as hpv
    T = sc.timer()
    print(f'hpvsim {hpv.__version__}')
    print(f'grid rel_sev_gt200={GRID}  n_agents={N_AGENTS} seeds={N_SEEDS} ms=5')
    print(f'baseline RWANDA_HIV_EFFECTS: {rs.rc.RWANDA_HIV_EFFECTS}')
    # HIV+ agg uses the fem-25+ denominator (apples-to-apples with the registry
    # adult-women target; HIV+ are all adults so it barely differs from all-fem).
    print(f'{"gt200":>7} {"HIV+ a35":>9} {"HIV+ a45":>9} {"HIV+ agg25":>11} '
          f'{"HIV- agg25":>11} {"cancers+":>9}')
    rows = []
    for gt in GRID:
        effects = _effects_with(gt)
        res_list = sc.parallelize(
            vf._run_one,
            iterkwargs=dict(seed=list(range(N_SEEDS))),
            kwargs=dict(n_agents=N_AGENTS, ms=5, effects=effects),
            serial=False,
        )
        rp, _, _ = vf._pool_byage(res_list, 'pos')
        hivp, npos, _ = vf._pool_agg(res_list, 'pos', 'nf_25')
        hivn, _, _ = vf._pool_agg(res_list, 'neg', 'nf_25')
        print(f'{gt:>7.2f} {rp[1]:>9.1f} {rp[2]:>9.1f} {hivp:>11.1f} '
              f'{hivn:>11.1f} {npos:>9.0f}', flush=True)
        rows.append((gt, hivp, hivn, rp.tolist()))

    best = min(rows, key=lambda r: abs(r[1] - TARGET_HIVP))
    print(f'\nclosest to HIV+ agg(25+) target {TARGET_HIVP}: '
          f'rel_sev_gt200={best[0]:.2f} -> HIV+ {best[1]:.1f}, HIV- {best[2]:.1f}')
    print(f'  HIV+ by-age at that point: {[round(x,1) for x in best[3]]} '
          f'(target {vf.TGT_BYAGE["pos"]})')
    T.toc('refit done')
