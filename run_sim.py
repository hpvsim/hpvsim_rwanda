"""
Define the HPVsim simulation for Rwanda
"""
#Set the working dir
import os
# os.chdir('D:\TxV_modeling\hpvsim_rwanda-main')

# Standard imports
import pylab as pl
import numpy as np
import sciris as sc
import hpvsim as hpv
import pandas as pd


from interventions import make_st, make_vx

# %% Settings and filepaths
# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True

# %% Simulation creation functions #calib=False: do not run default calibration
# calib=True: run default calibration
# calib_pars=None: use default calibration parameters
# calib_pars=calib_pars: use custom calibration parameters
def make_sim(calib=False, calib_pars=None, use_calib=True, debug=debug, add_vax=True, add_st=True, interventions=None,
            analyzers=None, seed=1, end=2100, hiv_pars=None):
    """
    Define the simulation
    """
    if end is None: end = 2100
    if calib: end = 2025
        
    # Basic parameters
    pars = sc.objdict(
        n_agents=[10e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=end,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location='rwanda',
        init_hpv_dist=dict(hpv16=0.4, hpv18=0.25, hi5=0.25, ohr=.1),
        init_hpv_prev={
            'age_brackets': np.array([12, 17, 24, 34, 44, 64, 80, 150]),     
            'm': np.array([0.0, 0.25, 0.6, 0.25, 0.05, 0.01, 0.0005, 0]),
            'f': np.array([0.0, 0.35, 0.7, 0.25, 0.05, 0.01, 0.0005, 0]),
        },
        ms_agent_ratio=100,                 #for every 1 "real" cancer case, 100 "normal" agents are simulated
        verbose=0.0,                        # the model runs silently without printing messages
        rand_seed=seed,
        model_hiv=True,
        hiv_pars=dict(),
    )

    # Sexual behavior parameters
    # Debut: derived by fitting to 2019-20 DHS
    # Women:
    # Age: 15,   18,   20,   22,   25
    # Prop_active: 2.1, 19.8, 41.7, 62.9, 81.7
    # Men:
    # Age:  15,   18,   20,   22,   25
    # Prop_active: 2.7, 14.4, 30.6, 47.6, 69.8
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    pars.debut = dict(
        f=dict(dist='lognormal', par1=20.96, par2=3.34),
        m=dict(dist='lognormal', par1=17.91, par2=2.83),
    )

    # Participation in marital and casual relationships
    # Derived to fit 2019-2020 DHS data
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    pars.layer_probs = dict(
        m=np.array([
            # Share of people of each age who are married
            [0, 5,    10,       15,      20,     25,      30,     35,      40,     45,    50,   55,   60,   65,   70,   75],
            [0, 0, 0.025,   0.0115,  0.1555,  0.313,  0.3875,  0.408,  0.3825,  0.334, 0.275, 0.20, 0.20, 0.20, 0.20, 0.20],  # Females
            [0, 0,  0.01,    0.023,   0.311,  0.626,   0.775,  0.816,   0.765,  0.668,  0.70, 0.80, 0.70, 0.60, 0.50, 0.60]]  # Males
        ),
        c=np.array([
            # Share of people of each age in casual partnerships
            [0, 5,  10,  15,  20,  25,  30,  35,  40,   45,   50,   55,   60,   65,   70,   75],
            [0, 0, 0.1, 0.6, 0.3, 0.2, 0.2, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0, 0, 0.1, 0.3, 0.4, 0.3, 0.3, 0.4, 0.5, 0.50, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]]
        ),
    )

    pars.m_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )
    pars.f_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )

    # HIV parameters and data
    pars.hiv_pars['art_failure_prob'] = 0.1
    hiv_datafile = ['data/hiv_incidence_rwanda.csv',
                    'data/rwanda_female_hiv_mortality.csv',
                    'data/rwanda_male_hiv_mortality.csv']
    art_datafile = ['data/rwanda_art_coverage_by_age_males.csv',
                    'data/rwanda_art_coverage_by_age_females.csv']

    # If calibration parameters have been supplied, use them here
    if calib_pars is None:
        # Use defaults
       if use_calib: calib_pars = sc.loadobj(f'results/rwanda_pars.obj')

    if hiv_pars is not None:
       pars.hiv_pars = sc.mergedicts(pars.hiv_pars, hiv_pars)

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    # Interventions
    interventions = sc.autolist(interventions)
    if add_vax: interventions += make_vx(end_year=end)
    if add_st: interventions += make_st(end_year=end)
            
    # Create the sim
    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=sc.tolist(analyzers),
                  hiv_datafile=hiv_datafile, art_datafile=art_datafile)

    return sim
                 
    
# %% Simulation running functions
def run_sim(
        analyzers=None, interventions=None, debug=debug, seed=1, verbose=1/4,
        do_save=False, end=None, add_vax=True, add_st=True, use_calib=True):

    # Make sim
    sim = make_sim(
        debug=debug,
        add_vax=add_vax,
        add_st=add_st,
        interventions=interventions,
        analyzers=analyzers,
        use_calib=use_calib,
        end=end,
    )
    sim.label = f'Sim--{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()        # Executes the simulation	Always — required to generate output
    sim.shrink()     # Minimizes memory by trimming details	Optional — for cleanup or saving

    if do_save:
        sim.save(f'results/rwanda.sim')

    return sim


def get_age_causal_df(sim=None):
    """
    Make age causal dataframe
    """
    dt_res = sim.get_analyzer('age_causal_infection')
    dt_dfs = sc.autolist()

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_causal)[np.array(dt_res.age_causal)<50]
    dt_df['Health event'] = 'Causal\ninfection'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cin)[np.array(dt_res.age_causal)<65]
    dt_df['Health event'] = 'HSIL'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cancer)[np.array(dt_res.age_causal)<90]
    dt_df['Health event'] = 'Cancer'
    dt_dfs += dt_df

    age_causal_df = pd.concat(dt_dfs)

    return age_causal_df


# %% Run as a script
if __name__ == '__main__':            # This is Python's standard way of saying: Only run the following code if this file is executed directly.But if another script imports this file (e.g., import my_script), the code in this block won’t run — which prevents unwanted side effects.

    T = sc.timer()

    # Make a list of what to run, comment out anything you don't want to run
    to_run = [
        'run_single',
    ]

    use_calib = True

    # Run and plot a single simulation
    # Takes <1min to run
    if 'run_single' in to_run:

        sim = run_sim(use_calib=use_calib, end=2035, analyzers=hpv.age_causal_infection(start_year=2020), debug=debug)  # Run the simulation
        # sim.plot()  # Plot the simulation
        # a = sim.get_analyzer()
        # a.plot()  # Save the plot to a file
        df = get_age_causal_df(sim)
        # save the dataframe to a file
        sc.saveobj('results/age_causal_infection.obj', df)


    T.toc('Done')  # Print out how long the run took



