from emcee.utils import MPIPool
import matplotlib.pyplot as plt
import subprocess as sp
import numpy as np
import pandas as pd
import seaborn as sns
import emcee

class MCMCrun:
    def __init__(self, name, nwalkers, burn_in=0):
        self.name = name 
        self.nwalkers = nwalkers
        self.burn_in = burn_in
        
        # read in chain from .csv
        self.chain = pd.read_csv(name + '/' + name + '_chain.csv')
    
        # indentify bad walkers and make clean chain
        last_step = self.chain.iloc[-nwalkers:]
        self.last_step = self.chain.iloc[-nwalkers:]
        self.bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % self.nwalkers
        self.clean_chain = self.chain.drop([row for row in self.chain.index 
            if row%self.nwalkers in self.bad_walkers])
        if len(self.bad_walkers !=0):
            print('walkers {} have not converged'.format(tuple(self.bad_walkers)))
        
        
    def evolution(self, show=False):
        plt.close()
        steps=self.chain.index // self.nwalkers

        color = 'red' if 0 in self.bad_walkers else 'black'
        axes = self.chain.iloc[0::self.nwalkers].plot(x=np.arange(steps[-1]+1),
            figsize=(7, 2.2*(len(self.chain.columns))), subplots=True, color=color, alpha=0.5)
            
        for i in range(self.nwalkers-1):
            color = 'red' if i+1 in self.bad_walkers else 'black'
            self.chain.iloc[i+1::self.nwalkers].plot(x=np.arange(steps[-1]+1),
                subplots=True, ax=axes, legend=False, color=color, alpha=0.5)

        # remove bad walkers before taking mean
        stepped_chain = self.clean_chain.set_index(steps)

        walker_means = pd.DataFrame([stepped_chain.loc[i].mean() for i in range(stepped_chain.index[-1])])
        walker_means.plot(subplots=True, ax=axes, legend=False, color='forestgreen', ls='--')

        plt.suptitle(self.name + ' walker evolution')
        plt.savefig('{}/evolution.png'.format(self.name), dpi=700)
        if show:
            plt.show()

# def corner_plot(run_name, nwalkers, stat_specs, burn_in=0, bad_walkers=[]):
#     """ Plot 'corner plot' of fit"""
#     plt.close()
# 
#     # read in samples
#     samples = pd.read_csv(run_name + '.csv')
#     rename(samples)
# 
#     # cut out burn in and bad walkers
#     # posterior = samples.iloc[-100:, :-1]
#     posterior = samples.iloc[burn_in*nwalkers:, :-1]
# 
#     last_step = samples.iloc[-nwalkers:]
#     bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % nwalkers
#     posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers], inplace=True)
#     print('walkers {} removed from posterior.'.format(tuple(bad_walkers)))
# 
#     # make corner plot
#     corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
#     corner.map_diag(sns.kdeplot, cut=0)
#     corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=True)
#     corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
#     corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)
# 
#     # get best_fit and posterior statistics
#     stats = posterior.describe().drop(['count', 'min', 'max'])
#     stats.loc['best fit'] = samples.drop('lnprob', 1).loc[samples['lnprob'].idxmax()]
#     stats = stats.iloc[[-1]].append(stats.iloc[:-1])
#     print(stats.round(2).to_latex())
# 
#     table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
#     table_ax.axis('off')
#     left, bottom = stat_specs
#     pd.plotting.table(table_ax, stats.round(2), bbox=[left, bottom, 1-left, .12], edges='open', colLoc='right')
# 
#     # hide upper triangle, so that it's a conventional corner plot
#     for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
#         corner.axes[i, j].set_visible(False)
# 
#     for ax in corner.axes.flat:
#         ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
#         ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))
# 
#     if burn_in == 0:
#         title = run_name + '.corner_ungroomed'
#     else:
#         title = run_name + '.corner_groomed'
# 
#     plt.subplots_adjust(top=0.9)
#     corner.fig.suptitle(r'{} Parameters, {} Walkers, {} Steps $\to$ {} Samples'
#         .format(posterior.shape[1], nwalkers, posterior.shape[0]//nwalkers, posterior.shape[0], fontsize=25))
#     corner.savefig(title + '.png', dpi=700)

def run_emcee(run_name, nsteps, nwalkers, lnprob, to_vary):
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
        
    # initiate sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name, to_vary), pool=pool)

    try:
        chain = pd.read_csv(run_name + '/chain.csv')
        start_step = chain.index[-1] // nwalkers
        print('Resuming {} at step {}'.format(run_name, start_step))
        pos = np.array(chain.iloc[-nwalkers:, :-1])
        with open(run_name + '/chain.csv', 'a') as f:
            f.write('\n')
    except IOError:
        
        sp.call(['mkdir', run_name])
        sp.call(['mkdir', run_name + '/model_files'])
        
        print('Starting {}'.format(run_name))
        start_step = 0
        
        with open(run_name + '/chain.csv', 'w') as f:
            np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob\n'),), 
                delimiter=',', fmt='%s')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
            for i in range(nwalkers)] 
            
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(start_step + i))
        pos, chisum, blob = result
        with open(run_name + '/chain.csv', 'a') as f: 
            np.savetxt(f, [np.append(pos[i], chisum[i]) for i in range(nwalkers)], delimiter=',')

    pool.close()
