import emcee
import numpy as np
import subprocess as sp
import os
import pandas as pd
from collections import OrderedDict
from disk_model import debris_disk, raytrace
from astropy.io import fits

class Observation:
    def __init__(self, root, name, rms):
        self.name = name
        self.root = root
        self.uvf  = fits.open(self.root + self.name + '.uvf')

        self.rms = rms

        self.dec = self.uvf[3].data['DECEPO'][0]
        self.ra = self.uvf[3].data['RAEPO'][0]

    def clean(self, show=True):
        """
        Clean and image (if desired) a observation-specific model.
        Either model image or residuals may be chosen.
        """

        print('')
        print('================================================================================')
        print('                            clean({})                                   ').format(self.name)
        print('================================================================================')

        # Set observation-specific clean filepath; clear filepaths
        sp.call('rm -rf {}.{{mp,bm,cl,cm}}'.format(self.root + self.name), shell=True)

        #Dirty clean; save rms for clean cutoff
        sp.call(['invert',
            'vis={}.vis'.format(self.root + self.name),
            'map={}.mp'.format(self.root + self.name),
            'beam={}.bm'.format(self.root + self.name),
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        imstat_out=sp.check_output(['imstat',
            'in={}.mp'.format(self.root + self.name),
            "region='boxes(256,0,512,200)'"])
        dirty_rms = float(imstat_out[-38:-29])
        print("Dirty rms is {}".format(dirty_rms))


        # Clean down to half the rms
        sp.call(['clean',
            'map={}.mp'.format(self.root + self.name),
            'beam={}.bm'.format(self.root + self.name),
            'out={}.cl'.format(self.root + self.name),
            'niters=100000', 'cutoff={}'.format(dirty_rms/2)])
        sp.call(['restor',
            'map={}.mp'.format(self.root + self.name),
            'beam={}.bm'.format(self.root + self.name),
            'model={}.cl'.format(self.root + self.name),
            'out={}.cm'.format(self.root + self.name)])

        # Display clean image with 2,4,6 sigma contours, if desired
        if show == True:

            # Display an unimportant imaage to get around the fact that the first
            # image displayed with cgdisp in a session can't be deleted
            sp.call(['cgdisp', 'in=cgdisp_start.im', 'type=p', 'device=/xs'])

            #Get rms for countours
            imstat_out = sp.check_output(['imstat',
                'in={}.cm'.format(self.root + self.name),
                "region='boxes(256,0,512,200)'"])
            clean_rms = float(imstat_out[-38:-29])
            print("Clean rms is {}".format(clean_rms))

            # Display
            sp.call(['cgdisp',
                'in={}.cm,{}.cm'.format(self.root + self.name, self.root + self.name),
                'type=p,c', 'device=/xs',
                'slev=a,{}'.format(clean_rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])



class Model:

    def make_fits(self, params):
        disk_params = params[:-1]
        PA = params[-1]

        model_disk = debris_disk.Disk(disk_params, obs=[300, 131, 300, 20])
        raytrace.total_model(model_disk,
            distance=9.91, # pc
            imres=0.03, # arcsec/pix
            xnpix=512, #image size in pixels
            freq0=self.observations[0][0].uvf[0].header['CRVAL4']*1e-9, # obs frequeency
            PA=PA,
            offs=[0.0,0.0], # offset from image center
            nchans=1, # continum
            isgas=False, # continuum!
            includeDust=True, #continuuum!!
            extra=0.0, # ?
            modfile = self.root + self.name)

    def obs_sample(self, obs, starflux):
        """
        Create model fits file with correct header information and sample using
        ALMA observation uv coverage to to create model .vis and .uvf files.
        """

        # define observation-specific model name and delete any preexisting conditions
        filepath = self.root + self.name + obs.name
        sp.call('rm -rf {}*'.format(filepath), shell=True)

        self.fits = fits.open(self.root + self.name + '.fits')

        # add starflux to central pixel
        crpix = int(self.fits[0].header['CRPIX1'])
        self.im = self.fits[0].data[0]
        self.im[crpix, crpix] += starflux

        # align with observation
        self.fits[0].header['CRVAL1'] = obs.ra
        self.fits[0].header['CRVAL2'] = obs.dec

        self.fits.writeto(filepath + '.fits')

        # Convert model into MIRIAD .im image file
        sp.call(['fits', 'op=xyin',
            'in={}.fits'.format(filepath),
            'out={}.im'.format(filepath)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

        # Sample the model image using the observation uv coverage
        sp.call(['uvmodel', 'options=replace',
            'vis={}.vis'.format(obs.root + obs.name),
            'model={}.im'.format(filepath),
            'out={}.vis'.format(filepath)], stdout=open(os.devnull, 'wb'))

        #Convert to UVfits
        sp.call(['fits', 'op=uvout',
            'in={}.vis'.format(filepath),
            'out={}.uvf'.format(filepath)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    def get_chi(self, obs):
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

        uvf  = fits.open(self.root + self.name + obs.name +'.uvf')
        modrlimwt = uvf[0].data['data']
        modrl_stokes = modrlimwt[::2, 0, 0, 0, 0, 0]
        modim_stokes = modrlimwt[::2, 0, 0, 0, 0, 1]

        # Calculate chi^2
        chi = np.sum((datrl_stokes - modrl_stokes)**2 * weights +
                     (datim_stokes - modim_stokes)**2 * weights)

        self.chis.append(chi)

    def clean(self, obs, residual=False, show=True):
        """
        Clean and image (if desired) a observation-specific model.
        Either model image or residuals may be chosen.
        """

        print('')
        print('================================================================================')
        if residual:
            print('                              clean({}.residual)                               ').format(obs.name)
        else:
            print('                                   clean({})                                   ').format(obs.name)
        print('================================================================================')

        # Set observation-specific clean filepath; clear filepaths
        filepath = self.root + self.name + '_' + obs.name
        if residual == True:
            filepath += '.residual'
        sp.call('rm -rf {}.{{mp,bm,cl,cm}}'.format(filepath), shell=True)

        # Clean down to half the observation rms
        sp.call(['invert',
            'vis={}.vis'.format(filepath),
            'map={}.mp'.format(filepath),
            'beam={}.bm'.format(filepath),
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        sp.call(['clean',
            'map={}.mp'.format(filepath),
            'beam={}.bm'.format(filepath),
            'out={}.cl'.format(filepath),
            'niters=100000', 'cutoff={}'.format(obs.rms/2)])
        sp.call(['restor',
            'map={}.mp'.format(filepath),
            'beam={}.bm'.format(filepath),
            'model={}.cl'.format(filepath),
            'out={}.cm'.format(filepath)])

        # Display clean image with 2,4,6 sigma contours, if desired
        if show == True:

            # Display an unimportant image to get around the fact that the first
            # image displayed with cgdisp in a session can't be deleted
            sp.call(['cgdisp', 'in=cgdisp_start.im', 'type=p', 'device=/xs'])

            #Get rms for countours
            imstat_out = sp.check_output(['imstat',
                'in={}.cm'.format(filepath),
                "region='boxes(256,0,512,200)'"])
            clean_rms = float(imstat_out[-38:-29])
            print("Clean rms is {}".format(clean_rms))

            # Display
            sp.call(['cgdisp',
                'in={}.cm,{}.cm'.format(filepath, filepath),
                'type=p,c', 'device=/xs',
                'slev=a,{}'.format(clean_rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])

    def residuals(self, obs, show=True):

        """
        Create model residuals (data - model), and clean//display if desired
        """

        #Set observation-specific filepath
        filepath = self.root + self.name + '_' + obs.name

        # Subtract model visibilities from data; outfile is residual visibilities
        sp.call(['uvmodel', 'options=subtract',
            'model={}.im'.format(filepath),
            'vis={}.vis'.format(obs.root + obs.name),
            'out={}.residual.vis'.format(filepath)])

        if show == True:
            self.clean(obs, residual=True)

    def __init__(self, params, observations, root, name=''):

        # assign name and set of observations
        self.name = name
        self.root = root
        self.observations = observations


        #split parameters into disk params and star params
        self.disk_params = params[:-len(observations)]
        self.starfluxes  = params[-len(observations):]

        self.make_fits(self.disk_params)

#==============================================================================#
# Create observations, default parameter dict
#==============================================================================#
# mar0 = Observation('aumic_band6_mar_spw0_FINAL', rms=6.5e-05)
# mar1 = Observation('aumic_band6_mar_spw1_FINAL', rms=6.124e-05)
# mar2 = Observation('aumic_band6_mar_spw2_FINAL', rms=6.068e-05)
# mar3 = Observation('aumic_band6_mar_spw3_FINAL', rms=6.468e-05)
# aug0 = Observation('aumic_band6_aug_spw0_FINAL', rms=5.879e-05)
# aug1 = Observation('aumic_band6_aug_spw1_FINAL', rms=5.336e-05)
# aug2 = Observation('aumic_band6_aug_spw2_FINAL', rms=6.092e-05)
# aug3 = Observation('aumic_band6_aug_spw3_FINAL', rms=5.558e-05)
# jun0 = Observation('aumic_band6_jun_spw0_FINAL', rms=5.369e-05)
# jun1 = Observation('aumic_band6_jun_spw1_FINAL', rms=4.658e-05)
# jun2 = Observation('aumic_band6_jun_spw2_FINAL', rms=5.083e-05)
# jun3 = Observation('aumic_band6_jun_spw3_FINAL', rms=5.559e-05)
# band6_observations=[[mar0, mar1, mar2, mar3],
#                     [aug0, aug1, aug2, aug3],
#                     [jun0, jun1, jun2, jun3]]
#
# params = OrderedDict([
#     ('temp_index',        -0.5),
#     ('m_disk',            3.67e-08),
#     ('sb_law',            2.3),
#     ('r_in',              8.8),
#     ('d_r',                31.5),
#     ('r_crit',            150.0),
#     ('inc',               89.5),
#     ('m_star',            0.31),
#     ('co_frac',           0.0001),
#     ('v_turb',            0.081),
#     ('Zq',                70.0),
#     ('column_densities', [0.79, 1000]),
#     ('abundance_bounds', [50, 500]),
#     ('hand',              -1),
#     ('rgrid_size',        500),
#     ('zgrid_size',        500),
#     ('l_star',            0.09),
#     ('scale_factor',      0.1),
#     ('pa',                128.41),
#     ('mar_starflux',      3.67e-4),
#     ('aug_starflux',      1.23e-4),
#     ('jun_starflux',      2.62e-4)])
