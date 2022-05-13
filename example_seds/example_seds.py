import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.coordinates import SkyCoord

from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

cm = 1/2.54

# TOPCSV='GLEAM-X J065228.65-255000.34.csv'
TOPCSV = 'GLEAM-X J055252.76-222513.62.csv'

markers = {
    'gleamx':{'marker':'.', 'color':'black', 'label':'GLEAM-X', 'markeredgewidth':0.1, 'elinewidth':0.5},
    'sumss':{'marker':'s', 'color':'green', 'label':'SUMSS', 'markersize':2},
    'nvss':{'marker':'+', 'color':'red', 'label':'NVSS'}

}

ref_nu = 200

class CoordStr:
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0:
            self.pos = args[0]
        else:
            self.pos = SkyCoord(*args, **kwargs)
        
        
    def _make_str(self, tex=False):
        _format = None if tex is False else 'latex'
        
        return (f'GLEAM-X '\
                f'J{self.pos.ra.to_string(unit=u.hourangle, sep="", precision=1, pad=True, format=_format)}' \
                f'{self.pos.dec.to_string(sep="", precision=0, alwayssign=True, pad=True, format=_format)}')

    def __str__(self):
        return self._make_str()
    
    def __repr__(self):
        return str(self)
     
    @property
    def tex(self):
        return self._make_str(tex=True)



def pl(nu, norm, alpha):
    return norm * nu **alpha 

def two_comp_model(nu, s_nu, alpha_thin, alpha_thick, nu_peak):
    spec_nu = nu / nu_peak
    
    c1 = s_nu / (1 - np.exp(-1))
    c2 = 1 - np.exp(-spec_nu**(alpha_thin - alpha_thick))
    c3 = spec_nu ** alpha_thick
          
    return c1 * c2 * c3 

def curved_law(nu, s_nu, alpha, q):
    spec_nu = nu / ref_nu
        
    return s_nu * spec_nu ** alpha * \
            np.exp(q * np.log(spec_nu)**2)

def fn_to_tex(fn, old=False):
    if old:
        # because I am a bad person
        return fn.replace('.csv','').replace('-','--')
    
    print(fn)
    print(fn[8:-4])
    sky_coord = SkyCoord(fn[8:-4], unit=(u.hourangle, u.degree))
    coord_str = CoordStr(sky_coord)

    return coord_str.tex


def overlay_box(ax, text, x=0.02, y=0.125):
    print(text)
    ax.text(
        x, 
        y,
        text, 
        transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.25')
    )

def make_ax1(ax1, nu, csv=None):
    if csv is None:
        ax1.plot(
            nu,
            pl(nu, 100., -0.8),
            lw=0.5,
        )
        title = 'GLEAM-XJ12345-234565'
    else:
        title = fn_to_tex(csv)
        df = pd.read_csv(csv)
        
        gx_mask = df.freq < 300

        ax1.errorbar(
            df.freq[gx_mask],
            df.flux[gx_mask],
            yerr=df.fluxerr[gx_mask],
            linestyle='None',
            **markers['gleamx']
        )

        sumss_mask = df.freq == 887.5
        ax1.errorbar(
            df.freq[sumss_mask],
            df.flux[sumss_mask],
            yerr=df.fluxerr[sumss_mask],
            linestyle='None',
            **markers['sumss']
        )

        nvss_mask = df.freq == 1400
        ax1.errorbar(
            df.freq[nvss_mask],
            df.flux[nvss_mask],
            yerr=df.fluxerr[nvss_mask],
            linestyle='None',
            **markers['nvss']
        )

        fit_p0 = [np.median(df.flux), -0.7, 1.5, 250]
        fit_res = curve_fit(
            two_comp_model,
            df.freq,
            df.flux,
            fit_p0,
            sigma=df.fluxerr,
            absolute_sigma=True
        )

        no_samps = 1000
        samps = np.random.multivariate_normal(
            fit_res[0], fit_res[1], size=no_samps
        ).swapaxes(0,1)
        
        freq = nu
        models = two_comp_model(
            nu[:, None],
            *samps
        )

        q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)
        
        ax1.plot(
            nu,
            q50,
            lw=0.5,
        )
        ax1.fill_between(
            nu,
            q16, q84, alpha=0.3
        )


    ax1.grid(
        which='both',
    )
    ax1.legend(loc='lower right')

    ax1.set(
        xscale='log',
        yscale='log',
        xlabel='Frequency (MHz)',
        ylabel='Flux density (mJy)',
    )
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))
#    ax1.xaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))

    overlay_box(
        ax1,
        f"{title}\nDouble Power Law",
        y=0.875,
        x=0.012
    )

def make_small_ax(ax, nu, xlabel=None, onright=False, csv=None, model=None):
    
    if csv is None:
        ax.plot(
            nu,
            pl(nu, 100., -0.8),
        )
        title = None
    else:
        title = fn_to_tex(csv)
        scale = 1000
        df = pd.read_csv(csv)
        ax.errorbar(
            df.freq,
            df.flux*scale,
            yerr=df.fluxerr*scale,
            linestyle='None',
            **markers['gleamx']
        )

    ax.grid(
        which='both',
    )

    ax.set(
        xscale='log',
        yscale='log',
        ylabel='Flux density (mJy)',
    )
    
    if model is not None:
        if model == 'pl':
            model = (pl, (np.median(df.flux), -0.7), 
            'Power Law')
        elif model == 'cpl':
            model = (curved_law, (np.median(df.flux), -0.7, 0), 'Curved Power Law')
        else:
            print(f"{model} is not currently supported")

        nu = np.geomspace(65, 235, 100)
        fit_func = model[0]
        fit_p0 = model[1]
        fit_res = curve_fit(
            fit_func,
            df.freq,
            df.flux*scale,
            fit_p0,
            sigma=df.fluxerr*scale,
            absolute_sigma=True
        )

        no_samps = 1000
        samps = np.random.multivariate_normal(
            fit_res[0], fit_res[1], size=no_samps
        ).swapaxes(0,1)
        
        freq = nu
        models = fit_func(
            nu[:, None],
            *samps
        )

        q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)
        
        ax.plot(
            nu,
            q50,
            lw=0.5,
        )
        ax.fill_between(
            nu,
            q16, q84, alpha=0.3
        )

    if onright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')

    
    if title is not None:
        overlay_box(
            ax,
            f"{title}\n{model[2]}",
            y=0.06,
            x=0.315 if onright is False else 0.02
        )

    if xlabel is None:
        print('--Setting to None')
        ax.set_xticklabels([])
        ax.set_xticks([])
    else:
        ax.set_xlabel(xlabel)    
    ax.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax.xaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax.yaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))



def make_sed_figure(output='test.png'):

    example_nu_large = np.geomspace(65, 1500, 200)
    example_nu = np.geomspace(65, 270, 200)

    ax1_loc = (0.1, 0.6, 0.8, 0.35)
    ax2_loc = (0.1, 0.315, 0.395, 0.2)
    ax3_loc = (0.505, 0.315, 0.395, 0.2)
    ax4_loc = (0.1, 0.075, 0.395, 0.23)
    ax5_loc = (0.505, 0.075, 0.395, 0.23)


    fig = plt.figure(figsize=(17*cm, 15*cm))

    ax1 = fig.add_axes(ax1_loc)
    make_ax1(ax1, example_nu_large, csv=TOPCSV)

    ax2 = fig.add_axes(ax2_loc)
    ax3 = fig.add_axes(ax3_loc)
    ax4 = fig.add_axes(ax4_loc)#, sharex=ax2)
    ax5 = fig.add_axes(ax5_loc)#, sharex=ax3)

    print('Top Left')
    make_small_ax(ax2, example_nu, csv='GLEAM-X J075203.35-211015.30.csv', model='cpl')

    print('Top Right')
    make_small_ax(ax3, example_nu, onright=True, csv='GLEAM-X J050107.00-304737.32.csv', model='pl')

    print('Bottom Left')
    make_small_ax(ax4, example_nu, xlabel='Frequency (MHz)', csv='GLEAM-X J134551.54-301504.30.csv', model='pl')

    print('Bottom Right')
    make_small_ax(ax5, example_nu, xlabel='Frequency (MHz)', onright=True, csv='GLEAM-X J052952.39-242742.42.csv', model='cpl')

    fig.savefig(output)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Make some nice SED figure')

    parser.add_argument('-o','--output', type=str, default='test.png', help='Name of the output file')

    args = parser.parse_args()

    make_sed_figure(
        output=args.output
    )

