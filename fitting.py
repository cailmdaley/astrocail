import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess as sp
from astropy.io import fits
import os
# from galario.double import chi2Image
from scipy.optimize import curve_fit

class Observation:
    def __init__(self, root, name,  rms):
        self.root = root
        self.name = name
        self.path = root + name

        self.rms = rms
        self.uvf  = fits.open(self.path + '.uvf')
        try:
            self.dec = self.uvf[0].data['OBSDEC'][0]
            self.ra = self.uvf[0].data['OBSRA'][0]
        except:
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
            #'fwhm=0.91',
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'], stdout=open(os.devnull, 'wb'))
       # imstat_out=sp.check_output(['imstat',
       #     'in={}.mp'.format(path),
       #     "region='boxes(256,0,512,200)'"])
       # dirty_rms = float(imstat_out[-38:-29])
       # print("Dirty rms is {} for {}".format(dirty_rms, path))

        sp.call(['clean',
            'map={}.mp'.format(path),
            'beam={}.bm'.format(path),
            'out={}.cl'.format(path),
            'niters=10000', 'cutoff={}'.format(rms/2.)], stdout=open(os.devnull, 'wb'))
        sp.call(['restor',
            'map={}.mp'.format(path),
            'beam={}.bm'.format(path),
            'model={}.cl'.format(path),
            'out={}.cm'.format(path)], stdout=open(os.devnull, 'wb'))
        imstat_out=sp.check_output(['imstat',
            'in={}.cm'.format(path),
            "region='boxes(256,0,512,200)'"])
        #rms = float(imstat_out[-38:-29])
        #print("Clean rms is {} for {}".format(rms, path))

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

    def get_chi(self, obs, suffix=''):
        """
        Return chi^2 statistics of model.
        """
        # array containing real/imaginary/weight values, squeeze empty dimensions
        data_rlimwt = obs.uvf[0].data['data'].squeeze()
        freq0 = obs.uvf[0].header['crval4']

        # get weights--same for both xx and yy, arbitrarily choose xx (index 0)
        weights = data_rlimwt[:,0,2]

        # real and imaginary arrays in Stokes I
        if data_rlimwt.shape[1] == 2: # polarized; turn to stokes
            data_real = (data_rlimwt[:,0,0]+data_rlimwt[:,1,0])/2.
            data_imag = (data_rlimwt[:,0,1]+data_rlimwt[:,1,1])/2.
        else: # already stokes
            data_real = data_rlimwt[:,0,0]
            data_imag = data_rlimwt[:,0,1]

        # open model uvf, get rlimwt array again
        model_uvf  = fits.open(self.path + suffix + '.uvf')
        model_rlimwt = (model_uvf[0].data['data']).squeeze()

        # get real and imaginary values, skipping repeating values created by uvmodel.
        # when uvmodel converts to Stokes I, it either puts the value in place
        # of BOTH xx and yy, or makes xx and yy the same.
        # either way, selecting every other value solves the problem.
        model_real = model_rlimwt[::2,0]
        model_imag = model_rlimwt[::2,1]

        # Calculate chi^2
        chi = np.sum((data_real - model_real)**2 * weights +
                     (data_imag - model_imag)**2 * weights)

        self.chis.append(chi)
        
    def get_chi_galario(self, obs, plot_resid=False):

        # - Read in observation visibilities
        freq0 = obs.uvf[0].header['crval4']
        uu = np.ascontiguousarray(obs.uvf[0].data['UU']*freq0, np.float64)
        vv = np.ascontiguousarray(obs.uvf[0].data['VV']*freq0, np.float64)
        vis_rlimwt = np.ascontiguousarray(obs.uvf[0].data['data'].squeeze(), np.float64)
        
        # real and imaginary arrays in Stokes I
        if vis_rlimwt.shape[1] == 2: # polarized; turn to stokes
            vis_real = (vis_rlimwt[:,0,0]+vis_rlimwt[:,1,0])/2.
            vis_imag = (vis_rlimwt[:,0,1]+vis_rlimwt[:,1,1])/2.
        else: # already stokes
            vis_real = vis_rlimwt[:,0,0]
            vis_imag = vis_rlimwt[:,0,1]
        # get weights--same for both xx and yy, arbitrarily choose xx (index 0)
        vis_weights = np.ascontiguousarray(vis_rlimwt[:,0,2])

        #Generate model visibilities
        model_fits  = fits.open(self.path + '.fits')
        model_image = np.ascontiguousarray(model_fits[0].data.squeeze(), np.float64)
        dxy = np.radians(np.abs(model_fits[0].header['cdelt1']))
        model_fits.close()
        
        chi = chi2Image(model_image, dxy, uu, vv, vis_real, vis_imag, vis_weights)
        self.chis.append(chi)

    def view_fits(self):
        model_image = fits.getdata(self.path + '.fits')[0]
        plt.imshow(model_image, origin='lower')
        plt.show(block=False)
    
        
        '''Calculate the raw chi-squared based on the difference between the model and data visibilities.
        :param datfile: (default = 'data/HD163296.CO32.regridded.cen15')
         The base name for the data file. The code reads in the visibilities from datfile+'.vis.fits'
         :param modfile" (default='model/testpy_alma')
         The base name for the model file. The code reads in the visibilities from modfile+'.model.vis.fits'
         :param new_weight:
         An array containing the weights to be used in the chi-squared calculation. This should have the same dimensions as the real and imaginary part of the visibilities (ie Nbas x Nchan)
         :param systematic:
         The systematic weight to be applied. The value sent with this keyword is used to scale the absolute flux level of the model. It is defined such that a value >1 decreases the model and a value <1 increases the model (the model visibilities are divided by the systematic parameter). It is meant to mimic a true flux in the data which is larger or smaller by the fraction systematic (e.g. specifying systematic=1.2 is equivalent to saying that the true flux of the data is 20% brighter than what has been observed, with this scaling applied to the model instead of changing the data)
         :param isgas:
         If the data is line emission then the data has an extra dimension covering the >1 channels. Set this keyword to ensure that the data is read in properly. Conversely, if you are comparing continuum data then set this keyword to False.
    '''
