from astropy.io import fits
from astropy.modeling import models, fitting
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, LinearLocator, AutoMinorLocator
from astrocail import colormaps
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set seaborn plot styles and color pallete
sns.set_style("ticks",
              {"xtick.direction": "in",
               "ytick.direction": "in"})
sns.set_context("talk")
cpal = colormaps.jesse_reds


class Figure:

    def get_fits(self, path):
        fits_file = fits.open(path)
        self.head = fits_file[0].header
        self.im = fits_file[0].data[0][0]
        # self.im[np.isnan(self.im)]=0.

        # change units to micro Jy
        self.im *= 1e6
        self.rms *= 1e6

        # Read in header spatial info to create ra
        nx = self.head['NAXIS1'];           ny = self.head['NAXIS2']
        xpix = self.head['CRPIX1'];         ypix = self.head['CRPIX2']
        xval = self.head['CRVAL1'];         yval = self.head['CRVAL2']
        self.xdelt = self.head['CDELT1'];   self.ydelt = self.head['CDELT2']
        
        # Convert from degrees to arcsecs
        self.ra_offset = np.array(((np.arange(nx) - xpix + 1) * self.xdelt) * 3600)
        self.dec_offset = np.array(((np.arange(ny) - ypix + 1) * self.ydelt) * 3600)


    def make_axis(self, ax):
        
        xmin = -5.0
        xmax = 5.0
        ymin = -5.0
        ymax = 5.0
        ax.set_xlim(xmax, xmin)
        ax.set_ylim(ymin, ymax)
        ax.grid(False)

        # Set x and y major and minor tics
        majorLocator = MultipleLocator(1)
        ax.xaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_locator(majorLocator)

        minorLocator = MultipleLocator(0.2)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)

        # Set x and y labels
        ax.set_xlabel(r'$\Delta \alpha$ (")', fontsize=18)
        ax.set_ylabel(r'$\Delta \delta$ (")', fontsize=18)
        ax.xaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=18)
        ax.yaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=18)
        ax.tick_params(which='both', right='on', labelsize=18)

        # Set labels depending on position in figure
        if np.where(self.axes == ax)[0] % self.width == 0: #left
            ax.tick_params(axis='y', labelright='off', right='on')
        elif np.where(self.axes == ax)[0] % self.width == self.width - 1: #right
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='y', labelleft='off', labelright='on')
        else: #middle
            ax.tick_params(axis='y', labelleft='off')
            ax.set_xlabel('')
            ax.set_ylabel('')

        # Set physical range of colour map
        self.extent = [self.ra_offset[0], self.ra_offset[-1],
                      self.dec_offset[-1], self.dec_offset[0]]


    def fill_axis(self, ax, cbspace):
        # Plot image as a colour map
        cmap = ax.imshow(self.im,
            extent=self.extent,
            vmin=np.min(self.im),
            vmax=np.max(self.im),
            cmap=cpal)

        
        if self.rms:
            # Set contour levels
            cont_levs = np.arange(3, 100, 3) * self.rms
            # add residual contours if resdiual exists; otherwise, add image contours
            try:
                ax.contour(self.resid,
                               levels=cont_levs,
                               colors='k',
                               linewidths=0.75,
                               linestyles='solid')
                ax.contour(self.resid,
                               levels=-1 * np.flip(cont_levs, axis=0),
                               colors='k',
                               linewidths=0.75,
                               linestyles='dashed')
            except AttributeError:
                ax.contour(self.ra_offset, self.dec_offset, self.im,
                               colors='k',
                               levels=cont_levs,
                               linewidths=0.75,
                               linestyles='solid')
                ax.contour(self.ra_offset, self.dec_offset, self.im,
                               levels=-1 * np.flip(cont_levs, axis=0),
                               colors='k',
                               linewidths=0.75,
                               linestyles='dashed')
        
        # Create the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="8%", pad=0.0)
        cbar = self.fig.colorbar(cmap, ax=ax, cax=cax, orientation='horizontal')
        cbar.ax.xaxis.set_tick_params(direction='out', length=3, which='major',
            bottom='off', top='on', labelsize=8, pad=-2, 
            labeltop='on', labelbottom='off')

        cbar.ax.xaxis.set_tick_params(direction='out', length=2, which='minor',
            bottom='off', top='on')
        
        tickmaj, tickmin = cbspace
        minorLocator = AutoMinorLocator(tickmaj / tickmin)
        cbar.ax.xaxis.set_minor_locator(minorLocator)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                                rotation=45, fontsize=18)
        cbar.set_ticks(np.arange(-10*tickmaj, 10*tickmaj, tickmaj))

        # Colorbar label
        cbar.ax.text(0.425, 0.320, r'$\muJy / beam$', fontsize=12,
                     path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

        # Overplot the beam ellipse
        try:
            beam_ellipse_color = 'k'
            bmin = self.head['bmin'] * 3600.
            bmaj = self.head['bmaj'] * 3600.
            bpa = self.head['bpa']

            el = Ellipse(xy=[4.2, -4.2], width=bmin, height=bmaj, angle=-bpa,
                edgecolor='k', hatch='///', facecolor='none', zorder=10)
            ax.add_artist(el)
        except KeyError: 
            pass

        # Plot the scale bar
        x = -3.015
        y = -4.7
        ax.plot(
            [x, x - 1],
            [y, y],
            '-', linewidth=2, color='k')
        ax.text(
            x + 0.32, y + 0.15, "10 au",
            fontsize=18,
            path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

        # Plot a cross at the source position
        # ax.plot([0.0], [0.0], '*', markersize=6, markeredgewidth=1, color='k')

        # # Add figure text
        # try:
        #     for t in self.text:
        #         ax.text(*t, fontsize=18,
        #             path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
        #             
        # except AttributeError:
        #     pass
            
    def quickview(self):
            plt.imshow(self.im, origin='lower')
            plt.show(block=False)

    def __init__(self, paths, rmses=None, cbspaces=(100.,20.)):
        
        self.paths = [paths] if type(paths) is not list else paths
        rmses = [rmses] if type(rmses) is not list else rmses
        # cbspaces = [cbspaces] if type(cbspaces) is not list else cbspaces
        
        try:
            self.height, self.width = np.shape(self.paths)
        except ValueError:
            self.width = np.shape(self.paths); self.height=1
        
        self.fig, self.axes = plt.subplots(self.height, self.width,     
            figsize=(11.6/2 * self.width, 6.5*self.height),
            sharex=False, sharey=False)
        plt.subplots_adjust(wspace=-0.0)
        
        for row in zip(self.axes, self.paths, rmses):
            for ax, path, rms in zip(*row):
                self.rms = rms
                self.get_fits(path)
                self.make_axis(ax)
                self.fill_axis(ax, cbspaces)
        plt.savefig('run7/figure.pdf', dpi=700)
            
            
            
            
# attempt at making kde plot
#=========================================================================
# run_name = 'run5_26walkers_10params'
# posterior = pd.read_csv(run_name + '.csv')
# 
# 
# multi = posterior[['m_disk', 'sb_law']]
# 
# sns.kdeplot(multi)
# plt.show()
# 
# bw_search = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0, multi.std().mean()/10, 5)}, cv=20)
# bw_search.fit(multi)
# multi.shape[0]**(-1./(multi.shape[1]+4))
# multi.std()
# **()
# 
# xx, yy = np.meshgrid(np.linspace(*multi.iloc[:,0].quantile([0,1]), num=100), 
#                      np.linspace(*multi.iloc[:,1].quantile([0,1]), num=100))
# test = np.array([xx.ravel(), yy.ravel()]).T
# kde=KernelDensity(bandwidth=multi.std().min()/10)
# kde.fit(multi)
# pdf = np.exp(kde.score_samples(multi)).reshape(len(xx), -1)
# 
# plt.contour(xx, yy, pdf )
# plt.savefig('test.png')
# plt.show()
# 
# 
# 
# def my_kde(df):
#     mdisk = posterior['sb_law'].values.reshape(-1,1)
#     bw_search = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0, mdisk.std(), 5)})
#     bw_search.fit(mdisk)
#     kde = bw_search.best_estimator_.fit((mdisk))
#     x = np.linspace(mdisk.min(), mdisk.max(), 1000)[:,None]
#     pdf = np.exp(kde.score_samples(x))
#     
#     plt.plot(x, pdf)
#     plt.show()
