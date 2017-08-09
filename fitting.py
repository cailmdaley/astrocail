from astropy.io import fits
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess as sp
import os
from emcee.utils import MPIPool
import matplotlib.pyplot as plt
import emcee

class Observation:
    def __init__(self, root, name,  rms):
        self.root = root
        self.name = name
        self.path = root + name
        
        self.rms = rms
        self.uvf  = fits.open(self.path + '.uvf')
        self.dec = self.uvf[3].data['DECEPO'][0]
        self.ra = self.uvf[3].data['RAEPO'][0]
        
    def clean(self, show=True):
        """
        Clean and image (if desired) a observation-specific model.
        Either model image or residuals may be chosen.
        """

        # Set observation-specific clean filepath; clear filepaths
        sp.call('rm -rf {}.{{mp,bm,cl,cm}}'.format(self.path), shell=True)

        #Dirty clean; save rms for clean cutoff
        sp.call(['invert',
            'vis={}.vis'.format(self.path),
            'map={}.mp'.format(self.path),
            'beam={}.bm'.format(self.path),
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        imstat_out=sp.check_output(['imstat',
            'in={}.mp'.format(self.path),
            "region='boxes(256,0,512,200)'"])
        dirty_rms = float(imstat_out[-38:-29])
        print("Dirty rms is {}".format(dirty_rms))


        # Clean down to half the rms
        sp.call(['clean',
            'map={}.mp'.format(self.path),
            'beam={}.bm'.format(self.path),
            'out={}.cl'.format(self.path),
            'niters=100000', 'cutoff={}'.format(dirty_rms/2)])
        sp.call(['restor',
            'map={}.mp'.format(self.path),
            'beam={}.bm'.format(self.path),
            'model={}.cl'.format(self.path),
            'out={}.cm'.format(self.path)])

        # Display clean image with 2,4,6 sigma contours, if desired
        if show == True:

            # Display an unimportant imaage to get around the fact that the first
            # image displayed with cgdisp in a session can't be deleted
            sp.call(['cgdisp', 'in=cgdisp_start.im', 'type=p', 'device=/xs'])

            #Get rms for countours
            imstat_out = sp.check_output(['imstat',
                'in={}.cm'.format(self.path),
                "region='boxes(256,0,512,200)'"])
            clean_rms = float(imstat_out[-38:-29])
            print("Clean rms is {}".format(clean_rms))

            # Display
            sp.call(['cgdisp',
                'in={}.cm,{}.cm'.format(self.path, self.path),
                'type=p,c', 'device=/xs',
                'slev=a,{}'.format(clean_rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])

class Model:

    def __init__(self, observations, root, name=''):
        self.observations = observations
        self.name = name
        self.root = root
        self.path = root + name
        self.chis = []
        
        self.delete() # delete any preexisting files that will conflict
        
    def delete(self):
        sp.call('rm -rf {}*'.format(self.path), shell=True)

        
    def obs_sample(self, obs, suffix=''):
        """
        Create model fits file with correct header information and sample using
        ALMA observation uv coverage to to create model .vis and .uvf files.
        """

        # define observation-specific model name and delete any preexisting conditions
        sp.call('rm -rf {}{{.im,.vis,.uvf}}'.format(self.path + suffix), shell=True)

        # Convert model into MIRIAD .im image file
        sp.call(['fits', 'op=xyin',
            'in={}.fits'.format(self.path),
            'out={}.im'.format(self.path + suffix)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

        # Sample the model image using the observation uv coverage
        sp.call(['uvmodel', 'options=replace',
            'vis={}.vis'.format(obs.path),
            'model={}.im'.format(self.path + suffix),
            'out={}.vis'.format(self.path + suffix)], stdout=open(os.devnull, 'wb'))

        #Convert to UVfits
        sp.call(['fits', 'op=uvout',
            'in={}.vis'.format(self.path + suffix),
            'out={}.uvf'.format(self.path + suffix)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            
    def make_residuals(self, obs, suffix='', show=False):

        """
        Create model residuals (data - model), and clean//display if desired
        """
        sp.call('rm -rf *.residuals.vis'.format(self.path + suffix), shell=True)

        # Subtract model visibilities from data; outfile is residual visibilities
        sp.call(['uvmodel', 'options=subtract',
            'vis={}.vis'.format(obs.path),
            'model={}.im'.format(self.path + suffix),
            'out={}.residuals.vis'.format(self.path + suffix)], stdout=open(os.devnull, 'wb'))
        
        if show == True:
            self.clean(obs, residual=True)

    def get_chi(self, obs, suffix=''):
        """
        Return chi^2 statistics of model.
        """

        datrlimwt = obs.uvf[0].data['data']
        datrl_xx = datrlimwt[:, 0, 0, 0, 0, 0, 0]
        datrl_yy = datrlimwt[:, 0, 0, 0, 0, 1, 0]
        datim_xx = datrlimwt[:, 0, 0, 0, 0, 0, 1]
        datim_yy = datrlimwt[:, 0, 0, 0, 0, 1, 1]
        weights =  datrlimwt[:, 0, 0, 0, 0, 0, 2]
        datrl_stokes = np.array((datrl_xx + datrl_yy) / 2.)
        datim_stokes = np.array((datim_xx + datim_yy) / 2.)

        uvf  = fits.open(self.path + suffix + '.uvf')
        modrlimwt = uvf[0].data['data']
        modrl_stokes = modrlimwt[::2, 0, 0, 0, 0, 0]
        modim_stokes = modrlimwt[::2, 0, 0, 0, 0, 1]

        # Calculate chi^2
        chi = np.sum((datrl_stokes - modrl_stokes)**2 * weights +
                     (datim_stokes - modim_stokes)**2 * weights)

        self.chis.append(chi)

    def clean(self, path, rms, show=True):
        """
        Clean and image (if desired) a miriad .vis
        """

        # Set observation-specific clean filepath; clear filepaths
        sp.call('rm -rf {}.{{mp,bm,cl,cm}}'.format(path), shell=True)

        # Clean down to half the observation rms
        sp.call(['invert',
            'vis={}.vis'.format(path),
            'map={}.mp'.format(path),
            'beam={}.bm'.format(path),
            'sup=0', 'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs'])
        sp.call(['clean',
            'map={}.mp'.format(path),
            'beam={}.bm'.format(path),
            'out={}.cl'.format(path),
            'cutoff={}'.format(rms/2.),
            'region=arcsec,box(-5,-5,5,5)', 'niters=10000'], stdout=open(os.devnull, 'wb'))
        sp.call(['restor',
            'map={}.mp'.format(path),
            'beam={}.bm'.format(path),
            'model={}.cl'.format(path),
            'out={}.cm'.format(path)])
            
        # Convert MIRIAD .im image file into fits
        sp.call(['fits', 'op=xyout',
            'in={}.cm'.format(path),
            'out={}.fits'.format(path)], stdout=open(os.devnull, 'wb'))
            
        # Display clean image with 2,4,6 sigma contours, if desired
        if show == True:

            # Display an unimportant image to get around the fact that the first
            # image displayed with cgdisp in a session can't be deleted
            sp.call(['cgdisp', 'in=cgdisp_start.im', 'type=p', 'device=/xs'])

            # Display
            sp.call(['cgdisp',
                'in={}.cm,{}.cm'.format(path, path),
                'type=p,c', 'device=/xs',
                'slev=a,{}'.format(rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])
            raw_input('\npress enter when ready to go on:')
            
    def view_fits(self):
        model_image = fits.getdata(self.path + '.fits')[0]
        plt.imshow(model_image, origin='lower')
        plt.show(block=False)



class MCMCrun:
    def __init__(self, name, nwalkers, burn_in=0):
        self.name = name 
        self.nwalkers = nwalkers
        self.burn_in = burn_in
        
        # read in chain from .csv
        self.chain = pd.read_csv(name + '/chain.csv')
    
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
