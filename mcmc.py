from emcee.utils import MPIPool
import matplotlib.pyplot as plt
import subprocess as sp
import numpy as np
import pandas as pd
import seaborn as sns
import emcee
from astrocail import plotting

class MCMCrun:
    def __init__(self, name, nwalkers, path=None, burn_in=0):
        
        # read in chain from .csv
        self.main = pd.read_csv(path + '.csv') if path else pd.read_csv(name + '/' + name + '_chain.csv')
        self.name = name 
        self.nwalkers = nwalkers
        self.burn_in = burn_in
        self.nsteps = self.main.shape[0] / nwalkers
        
        last_steps = self.main.iloc[self.nsteps/2:]
        last_steps.index %= nwalkers
        self.bad_walkers = []
        for i in range(nwalkers):
            if last_steps.loc[i].duplicated(keep=False).sum() == len(last_steps) / nwalkers:
                self.bad_walkers.append(i)
                
                
        self.burnt_in = self.main.iloc[burn_in*nwalkers:]
        self.converged = self.burnt_in.drop([row for row in self.burnt_in.index 
            if row%nwalkers in self.bad_walkers])
        print('Removed bad walkers {}.'.format(tuple(self.bad_walkers)))
        print('Removed burn-in phase (first {} steps).'.format(burn_in))
        
        self.groomed = self.converged[self.converged['lnprob'] != -np.inf]
        
    def evolution(self, show=False):
        print('Making walker evolution plot...')
        plt.close()
        
        main = self.main.copy()
        main['lnprob'][np.isinf(main['lnprob'])] = np.nan
        
        

        color = 'red' if 0 in self.bad_walkers else 'black'
        axes = main.iloc[0::self.nwalkers].plot(
            x=np.arange(self.nsteps), figsize=(7, 2.2*(len(main.columns))), 
            subplots=True, color=color, alpha=0.5)
            
        for i in range(self.nwalkers-1):
            color = 'red' if i+1 in self.bad_walkers else 'black'
            main.iloc[i+1::self.nwalkers].plot(
                x=np.arange(self.nsteps), subplots=True, ax=axes, 
                legend=False, color=color, alpha=0.5)

        
        main[np.isnan(main['lnprob'])] = np.nan
        print(main[main['lnprob'] == np.nan])
        main.index //= self.nwalkers
        
        walker_means = pd.DataFrame([main.loc[i].mean() for i in range(self.nsteps)])
        walker_means.plot(subplots=True, ax=axes, legend=False, color='forestgreen', ls='--')
        
        plt.suptitle(self.name + ' walker evolution')
        plt.savefig(self.name + '/' + self.name + '_evolution.pdf'.format(self.name), dpi=700)
        if show:
            plt.show()

    def kde(self, show=False):
        print('Generating posterior kde plots...')
        plt.close()
        
        nrows, ncols = (2, int(np.ceil(self.groomed.shape[1]/2.)))
        fig, axes = plt.subplots(nrows, ncols, figsize=(2.5*ncols, 2.5*nrows))
        
        # plot kde of each free parameter
        for i, param in enumerate(self.groomed.columns):
            ax = axes.flatten()[i]
            
            for tick in ax.get_xticklabels(): tick.set_rotation(30)    
            ax.set_title(param)
            ax.tick_params(axis='y', left='off', labelleft='off')
            
            samples = self.groomed[param]
            plotting.my_kde(samples, ax=ax)
            
            percentiles = samples.quantile([.16,.5,.84])
            ax.axvline(percentiles.iloc[0], lw=1, ls='dotted', color='k', alpha=0.5)
            ax.axvline(percentiles.iloc[1], lw=1.5, ls='dotted', color='k')
            ax.axvline(percentiles.iloc[2], lw=1, ls='dotted', color='k', alpha=0.5)
            
            
            
            
        # bivariate kde to fill last subplot
        ax = axes.flatten()[-1]
        for tick in ax.get_xticklabels(): tick.set_rotation(30)    
        sns.kdeplot(self.groomed[r'$i$ ($\degree$)'], self.groomed[r'Scale Factor'], shade=True, cmap='Blues', n_levels=6, ax=ax);
        ax.tick_params(axis='y', left='off', labelleft='off', right='on', labelright='on')
        
        
        # adjust spacing and save
        plt.tight_layout()
        plt.savefig(self.name + '/' + self.name + '_kde.pdf'.format(self.name), dpi=700)
        
        # if show:
        plt.show(block=False)
        
        
def corner(run_name, nwalkers, stat_specs, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()

    # make corner plot
    corner = sns.PairGrid(self.groomed, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=True)
    corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)

    # get best_fit and posterior statistics
    stats = self.groomed.describe().drop(['count', 'min', 'max'])
    stats.loc['best fit'] = samples.drop('lnprob', 1).loc[samples['lnprob'].idxmax()]
    stats = stats.iloc[[-1]].append(stats.iloc[:-1])
    # print(stats.round(2).to_latex())

    table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
    table_ax.axis('off')
    left, bottom = stat_specs
    pd.plotting.table(table_ax, stats.round(2), bbox=[left, bottom, 1-left, .12], edges='open', colLoc='right')

    # hide upper triangle, so that it's a conventional corner plot
    for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
        corner.axes[i, j].set_visible(False)

    for ax in corner.axes.flat:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))

        title = self.name + 'Corner Plot'

    plt.subplots_adjust(top=0.9)
    corner.fig.suptitle(r'{} Parameters, {} Walkers, {} Steps $\to$ {} Samples'
        .format(posterior.shape[1], nwalkers, posterior.shape[0]//nwalkers, posterior.shape[0], fontsize=25))
    plt.savefig(self.name + '/' + self.name + '_corner.pdf'.format(self.name), dpi=700)

def run_emcee(run_name, nsteps, nwalkers, lnprob, to_vary):
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
        
    # initiate sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name, to_vary), pool=pool)

    try:
        chain = pd.read_csv(run_name + '/' + run_name + '_chain.csv')
        start_step = chain.index[-1] // nwalkers
        print('Resuming {} at step {}'.format(run_name, start_step))
        pos = np.array(chain.iloc[-nwalkers:, :-1])
        with open(run_name + '/' + run_name + '_chain.csv', 'a') as f:
            f.write('\n')
    except IOError:
        
        sp.call(['mkdir', run_name])
        sp.call(['mkdir', run_name + '/model_files'])
        
        print('Starting {}'.format(run_name))
        start_step = 0
        
        with open(run_name + '/' + run_name + '_chain.csv', 'w') as f:
            np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob'),), 
                delimiter=',', fmt='%s')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
            for i in range(nwalkers)] 
            
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(start_step + i))
        pos, chisum, blob = result
        with open(run_name + '/' + run_name + '_chain.csv', 'a') as f: 
            np.savetxt(f, [np.append(pos[i], chisum[i]) for i in range(nwalkers)], delimiter=',')

    pool.close()
