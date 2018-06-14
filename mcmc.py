from emcee.utils import MPIPool
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import subprocess as sp
import sys
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_style('ticks')
import emcee
from astrocail import plotting
import time

class MCMCrun:
    def __init__(self, name, nwalkers, path=None, burn_in=0):

        # read in chain from .csv
        self.main = pd.read_csv(path + '.csv') if path else pd.read_csv(name + '/' + name + '_chain.csv')
        self.name = name
        self.nwalkers = nwalkers
        self.nsteps = self.main.shape[0] // nwalkers

        # groom the chain
        self.burnt_in = self.main.iloc[burn_in*nwalkers:, :]
        self.groomed = self.burnt_in.loc[self.burnt_in.loc[:,'lnprob'] != -np.inf, :]
        print('Removed burn-in phase (step 0 through {} ).'.format(burn_in))

    def evolution(self):
        print('Making walker evolution plot...')
        plt.close()

        stepmin, stepmax = 0, self.nsteps
    
        main = self.main.copy().iloc[stepmin * self.nwalkers:stepmax * self.nwalkers,:]

        axes = main.iloc[0::self.nwalkers].plot(
            x=np.arange(stepmin,stepmax), figsize=(7, 2.0*(len(main.columns))),
            subplots=True, color='black', alpha=0.1)

        for i in range(self.nwalkers-1):
            main.iloc[i+1::self.nwalkers].plot(
                x=np.arange(stepmin, stepmax), subplots=True, ax=axes,
                legend=False, color='black', alpha=0.1)

            # make y-limits on lnprob subplot reasonable
            axes[-1].set_ylim(main.iloc[-1 * self.nwalkers:, -1].min(), main.lnprob.max())

        # if you want mean at each step over plotted:
        # main.index //= self.nwalkers
        # walker_means = pd.DataFrame([main.loc[i].mean() for i in range(self.nsteps)])
        # walker_means.plot(subplots=True, ax=axes, legend=False, color='forestgreen', ls='--')

        plt.suptitle(self.name + ' walker evolution')
        plt.savefig(self.name + '/' + self.name + '_evolution.png'.format(self.name))#, dpi=1)

    def kde(self):
        print('Generating posterior kde plots...')
        plt.close()

        nrows, ncols = (2, int( np.ceil( (self.groomed.shape[1] - 1) / 2. ) ) )
        fig, axes = plt.subplots(nrows, ncols, figsize=(2.5*ncols, 2.5*nrows))

        # plot kde of each free parameter
        for i, param in enumerate(self.groomed.columns[:-1]):
            ax = axes.flatten()[i]

            for tick in ax.get_xticklabels(): tick.set_rotation(30)
            ax.set_title(param)
            ax.tick_params(axis='y', left='off', labelleft='off')

            samples = self.groomed.loc[:,param]
            plotting.my_kde(samples, ax=ax)

            percentiles = samples.quantile([.16,.5,.84])
            ax.axvline(percentiles.iloc[0], lw=1, ls='dotted', color='k', alpha=0.5)
            ax.axvline(percentiles.iloc[1], lw=1.5, ls='dotted', color='k')
            ax.axvline(percentiles.iloc[2], lw=1, ls='dotted', color='k', alpha=0.5)


        # hide unfilled axes
        for ax in axes.flatten()[self.groomed.shape[1]:]:
            ax.set_axis_off()

        # bivariate kde to fill last subplot
        # ax = axes.flatten()[-1]
        # for tick in ax.get_xticklabels(): tick.set_rotation(30)
        # sns.kdeplot(self.groomed[r'$i$ ($\degree$)'], self.groomed[r'Scale Factor'], shade=True, cmap='Blues', n_levels=6, ax=ax);
        # ax.tick_params(axis='y', left='off', labelleft='off', right='on', labelright='on')

        # adjust spacing and save
        plt.tight_layout()
        plt.savefig(self.name + '/' + self.name + '_kde.png'.format(self.name))
        plt.show()


    def corner(self, variables=None):
        """ Plot 'corner plot' of fit"""
        plt.close()

        # get best_fit and posterior statistics
        stats = self.groomed.describe(percentiles=[0.16,0.84]).drop(['count', 'min', 'max', 'mean'])
        stats.loc['best fit'] = self.main.loc[self.main['lnprob'].idxmax()]
        stats = stats.iloc[[-1]].append(stats.iloc[:-1])
        stats.loc[['16%','84%'], :] -= stats.loc['50%',:]
        stats = stats.reindex(['50%', '16%', '84%', 'best fit', 'std'], copy=False)
        print(stats.T.round(6).to_string())

        # make corner plot
        corner = sns.PairGrid(data=self.groomed, diag_sharey=False, despine=False,
            vars=variables)
            
        if variables is not None: 
            corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.1)
        else:
            corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=18, shade=True)
        
        corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)
        corner.map_diag(sns.kdeplot, cut=0)
        
        if variables is None:
            # get best_fit and posterior statistics
            stats = self.groomed.describe(percentiles=[0.16,0.84]).drop(['count', 'min', 'max', 'mean'])
            stats.loc['best fit'] = self.main.loc[self.main['lnprob'].idxmax()]
            stats = stats.iloc[[-1]].append(stats.iloc[:-1])
            stats.loc[['16%','84%'], :] -= stats.loc['50%',:]
            stats = stats.reindex(['50%', '16%', '84%', 'best fit', 'std'], copy=False)
            print(stats.T.round(3).to_string())
            # print(stats.round(2).to_latex())

            # add stats to corner plot as table
            table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
            table_ax.axis('off')
            left, bottom = 0.15, 0.83
            pd.plotting.table(table_ax, stats.round(2), bbox=[left, bottom, 1-left, .12], edges='open', colLoc='right')
            
            corner.fig.suptitle(r'{} Parameters, {} Walkers, {} Steps $\to$ {} Samples'
                .format(self.groomed.shape[1], self.nwalkers,
                self.groomed.shape[0]//self.nwalkers, self.groomed.shape[0],
                fontsize=25))
            tag=''
        else:
            tag = '_subset'

        # hide upper triangle, so that it's a conventional corner plot
        for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
            corner.axes[i, j].set_visible(False)

        # fix decimal representation
        for ax in corner.axes.flat:
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))

        plt.subplots_adjust(top=0.9)
        plt.savefig(self.name + '/' + self.name + tag + '_corner.png'.format(self.name))

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
        end = np.array(chain.iloc[-nwalkers:,:])
        print('Start step: {}'.format(np.mean(end[:,-1])))
        
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

    lnprobs = []
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        old_lnprobs = np.copy(lnprobs)
        pos, lnprobs, blob = result
        print("Step {}: {}".format(start_step + i, np.mean(lnprobs)))
        # print('Acceptances: {}'.format([lnprob for lnprob in lnprobs if lnprob not in old_lnprobs]))
        # print('')
        # print(lnprobs)
        # print(np.mean(pos))
        with open(run_name + '/' + run_name + '_chain.csv', 'a') as f:
            np.savetxt(f, [np.append(pos[k], lnprobs[k]) for k in range(nwalkers)], delimiter=',')
                

    pool.close()

def run_emcee_simple(run_name, nsteps, nwalkers, lnprob, to_vary, burn_in=0, pool=False, resume=False):
    if pool:
        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    
    start = time.time()
    # initiate sampler chain
    ndim = len(to_vary[0])
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)

    if resume:
        chain = pd.read_csv(run_name + '/' + run_name + '_chain.csv')
        start_step = chain.index[-1] // nwalkers
        print('Resuming {} at step {}'.format(run_name, start_step))
        
        with open(run_name + '/' + run_name + '_chain.csv', 'a') as f:
            f.write('\n')
        pos = np.array(chain.iloc[-nwalkers:, :-1])
        
    else:
        sp.call('rm -rf ' + run_name, shell=True)
        sp.call(['mkdir', run_name])
        print('Starting {}'.format(run_name))
        start_step = 0
    
        with open(run_name + '/' + run_name + '_chain.csv', 'w') as f:
            f.write(','.join([param[0] for param in to_vary] + ['lnprob']) + '\n')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
                for i in range(nwalkers)] 
                
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(start_step + i))
        pos, chisum, blob = result
        with open(run_name + '/' + run_name + '_chain.csv', 'a') as f:
            for i in range(nwalkers):
                f.write(','.join(map(str, np.append(pos[i], chisum[i]))) + '\n')
    print('{} samples in {:.1f} seconds'.format(nsteps*nwalkers, time.time() - start))   
    
    if pool:
        pool.close()
    
    return MCMCrun(run_name, nwalkers, burn_in=burn_in)
