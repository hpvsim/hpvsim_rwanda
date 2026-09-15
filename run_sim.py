"""
Define the HPVsim simulation for Rwanda
"""
# Standard imports
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv


from interventions import make_st, make_vx

# %% Settings and filepaths
# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True

# HIV/ART inputs, in the four fixed filenames hpv.data.load_hiv_data expects.
hiv_datafolder = 'data/hiv'


# %% Simulation creation functions #calib=False: do not run default calibration
# calib=True: run default calibration
# calib_pars=None: use default calibration parameters
# calib_pars=calib_pars: use custom calibration parameters
def make_sim(calib=False, calib_pars=None, use_calib=True, debug=debug, add_vax=True, add_st=True, interventions=None,
            analyzers=None, seed=1, stop=2100, hiv_pars=None):
    """
    Define the simulation
    """
    if stop is None: stop = 2100
    if calib: stop = 2025

    # Runtime pars
    pars = sc.objdict(
        ms_agent_ratio=100,
        verbose=0.0,
        debut_f=ss.normal(mean=21.0, std=3.34),
        debut_m=ss.normal(mean=22.5, std=2.83),
        layer_probs_marital=np.array([
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75],
            [0.0000, 0.0000, 0.0963, 0.0452, 0.4914, 0.7772, 0.8593, 0.8772, 0.8546, 0.8033, 0.7237, 0.5904, 0.5904, 0.5904, 0.5904, 0.5904],
            [0.0000, 0.0000, 0.0394, 0.0889, 0.7746, 0.9804, 0.9974, 0.9989, 0.9970, 0.9879, 0.9919, 0.9984, 0.9919, 0.9744, 0.9375, 0.9744]]),
        layer_probs_casual=np.array([
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75],
            [0.0000, 0.0000, 0.3439, 0.9744, 0.7599, 0.5904, 0.5904, 0.5904, 0.5904, 0.1855, 0.0394, 0.0394, 0.0394, 0.0394, 0.0394, 0.0394],
            [0.0000, 0.0000, 0.3439, 0.7599, 0.8704, 0.7599, 0.7599, 0.8704, 0.9375, 0.9375, 0.0394, 0.0394, 0.0394, 0.0394, 0.0394, 0.0394]]),
        m_partners_marital=ss.poisson(lam=0.01),
        m_partners_casual=ss.poisson(lam=0.2),
        f_partners_marital=ss.poisson(lam=0.01),
        f_partners_casual=ss.poisson(lam=0.2),
    )

    if calib_pars is None and use_calib:
        calib_pars = sc.loadobj('results/rwanda_pars.obj')
    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)
        # calib was run with reseed=True; the fitted seed lives in pars and
        # always wins so the calibrated fit is exactly reproducible.
        if 'rand_seed' in pars:
            seed = int(pars.pop('rand_seed'))

    # v2 used `art_failure_prob=0.1`; v3 complement is `p_effective_art`
    # (Bernoulli p that an ART recipient becomes virally suppressed).
    base_hiv = dict(p_effective_art=0.9)
    hiv_pars = sc.mergedicts(base_hiv, hiv_pars)

    interventions = sc.autolist(interventions)
    if add_vax: interventions += make_vx(end_year=stop)
    if add_st: interventions += make_st(end_year=stop)

    return hpv.Sim(
        pars=pars,
        location='rwanda',
        genotypes=[16, 18, 'hi5', 'ohr'],
        init_hpv_dist=dict(hpv16=0.4, hpv18=0.25, hi5=0.25, ohr=.1),
        n_agents=[10e3, 1e3][debug],
        start=[1960, 1980][debug],
        stop=stop,
        dt=[0.25, 1.0][debug],
        rand_seed=seed,
        interventions=interventions,
        analyzers=sc.tolist(analyzers),
        model_hiv=True,
        hiv_data=hiv_datafolder,
        hiv_pars=hiv_pars,
    )

    return sim


# %% Simulation running functions
def run_sim(
        analyzers=None, interventions=None, debug=debug, seed=1, verbose=1/4,
        do_save=False, stop=None, add_vax=True, add_st=True, use_calib=True):

    # Make sim
    sim = make_sim(
        debug=debug,
        add_vax=add_vax,
        add_st=add_st,
        interventions=interventions,
        analyzers=analyzers,
        use_calib=use_calib,
        stop=stop,
    )
    sim.label = f'Sim--{seed}'

    # Run
    sim.pars.verbose = verbose
    sim.run()        # Executes the simulation
    sim.shrink()     # Minimizes memory by trimming details

    if do_save:
        sim.save(f'results/rwanda.sim')

    return sim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    # Make a list of what to run, comment out anything you don't want to run
    to_run = [
        'run_single',
    ]

    use_calib = True

    # Run and plot a single simulation (<1 min)
    if 'run_single' in to_run:
        sim = run_sim(use_calib=use_calib, stop=2025, debug=debug)

    T.toc('Done')



