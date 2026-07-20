"""Quick diagnostic: female HPV prevalence + total cancer ASR trajectory,
to separate a network collapse (low HPV prevalence) from a natural-history
cancer-level regression (HPV prevalence OK but cancer low)."""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')
import numpy as np
import run_sim as rs

sim = rs.make_sim(add_vax=False, add_st=False, seed=0, n_agents=15000,
                  start=1960, stop=2020, dt=0.25, ms_agent_ratio=5)
sim.run()

# ASR from the RwandaReport analyzer
rep = next(a for a in sim.analyzers.values() if isinstance(a, rs.RwandaReport))
tab = rep.annual_table()
print('year  ASR   inc_hiv  inc_nohiv')
for y in [1990, 2000, 2010, 2015, 2019]:
    idx = np.where(tab['year'] == y)[0]
    if len(idx):
        i = int(idx[0])
        print(f'{y}: {tab["asr_cancer_incidence"][i]:6.1f} '
              f'{tab["cancer_incidence_with_hiv"][i]:8.1f} '
              f'{tab["cancer_incidence_no_hiv"][i]:8.1f}')

# Female HPV prevalence (any genotype) among adult women, scale-weighted, at end
from hpvsim.hpv import HPV
p = sim.people
w = p.scale.values
fem = p.female.values & p.alive.values
adult = fem & (p.age.values >= 15) & (p.age.values < 50)
anyhpv = np.zeros(len(w), bool)
for m in sim.diseases.values():
    if isinstance(m, HPV):
        anyhpv |= m.infected.values
prev = (w * (anyhpv & adult)).sum() / (w * adult).sum()
print(f'\nfemale 15-49 any-HPV prevalence at 2019: {prev*100:.1f}%')
