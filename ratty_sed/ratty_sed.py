import os
from argparse import ArgumentParser
import logging
from collections import Counter
from multiprocessing import Pool
import warnings 

import matplotlib
matplotlib.use('Agg')

from reproject.mosaicking import find_optimal_celestial_wcs
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
import astropy.units as u 
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats.circstats import circmean
from astropy.nddata import Cutout2D
from astropy.io import fits

warnings.simplefilter('ignore', category=AstropyWarning)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

cm = lambda a: a/2.54

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

OUTDIR='SEDs'
CHANSN = ['072_080MHz', '080_088MHz', '088_095MHz', '095_103MHz', '103_111MHz', '111_118MHz', '118_126MHz',
 '126_134MHz', '139_147MHz', '147_154MHz', '154_162MHz', '162_170MHz', '170_177MHz', 
 '177_185MHz', '185_193MHz', '193_200MHz', '200_208MHz', '208_216MHz', '216_223MHz', '223_231MHz']

GLEAMX_FREQS = np.array([np.mean([float(i) for i in name.replace('MHz','').split('_')[-2:]]) for name in CHANSN])
GLEAMX_INT = [f'int_flux_N_{c}' for c in CHANSN]
GLEAMX_ERR = [f'err_int_flux_N_{c}' for c in CHANSN]

class CoordStr:
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0:
            self.pos = args[0]
        else:
            self.pos = SkyCoord(*args, **kwargs)
        
        
    def _make_str(self, tex=False):
        _format = None if tex is False else 'latex'
        
        return (f'GLEAM-X '\
                f'J{self.pos.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True, format=_format)}' \
                f'{self.pos.dec.to_string(sep="", precision=2, alwayssign=True, pad=True, format=_format)}')

    def __str__(self):
        return self._make_str()
    
    def __repr__(self):
        return str(self)
     
    @property
    def tex(self):
        return self._make_str(tex=True)


def get_freq_flux_err(row, apply_mask = True, internal_scale=0.02):
    
    freq = GLEAMX_FREQS
    int_flux = np.array([row[i] for i in GLEAMX_INT])
    err_flux = np.array([row[i] for i in GLEAMX_ERR])
    
    # Gleam-x internal flux error
    err_flux = np.sqrt( err_flux ** 2 + (int_flux * internal_scale)**2)
        
    if apply_mask:
        mask = np.isfinite(int_flux) & (int_flux > 0) \
            & np.isfinite(err_flux) & (err_flux > 0)
        
        freq = freq[mask]
        int_flux = int_flux[mask]
        err_flux = err_flux[mask]
    
    return freq, int_flux*1000, err_flux*1000

def add_beam(ax, wcs, psf_file, psf_pos, loc, hatch):
    
    with fits.open(psf_file) as psf_fits:
        psf_wcs = WCS(psf_fits[0].header).celestial
        x_pix, y_pix = psf_wcs.all_world2pix(psf_pos.ra, psf_pos.dec, 0)
        beamsize = psf_fits[0].data[:3, int(y_pix), int(x_pix)]

    pix_scale = proj_plane_pixel_scales(wcs.celestial)
    
    sx = pix_scale[0]
    sy = pix_scale[1] 

    degrees_per_pixel = np.sqrt(sx * sy)
    
    bmaj = beamsize[0] / degrees_per_pixel
    bmin = beamsize[1] / degrees_per_pixel
    bpa = beamsize[2]

    bmaj = beamsize[0] / sx
    bmin = beamsize[1] / sy
    
    beam1 = AnchoredEllipse(ax.transData, width=bmaj, height=bmin, angle=bpa,
                         loc=loc, pad=0.15, borderpad=0.4,
                         frameon=False)
    
    beam1.ellipse.set_fc('white')
    beam1.ellipse.set_hatch('//')
    beam1.ellipse.set_fill('white')
    # beam1.ellipse.set_alpha(0.5)
    ax.add_artist(beam1)


def overlay_box(ax, text, x=0.02, y=0.125):
    ax.text(
        x, 
        y,
        text, 
        transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.25')
    )

def plot_img_sed(idx, isl_df, img_deep_path, img_low_path, sep, deep_psf, low_psf):
    if len(isl_df) != 2:
        return 

    comp_pos = SkyCoord(
        isl_df.ref_ra * u.deg, 
        isl_df.ref_dec * u.deg
    )

    if comp_pos[0].separation(comp_pos[1]) > sep:
        return

    mean_pos = SkyCoord(circmean(comp_pos.ra), np.mean(comp_pos.dec))
    img_fac = 15
    with fits.open(img_deep_path) as img_fits:
        w = WCS(img_fits[0].header)
        cutout = Cutout2D(
            img_fits[0].data,
            mean_pos,
            sep*img_fac,
            wcs=w,
            copy=True
        )
        new_wcs, new_shape = find_optimal_celestial_wcs(
            ((cutout.data, cutout.wcs),), 
            projection='SIN',
            reference=mean_pos
        )
        


    with fits.open(img_low_path) as img_fits:
        w = WCS(img_fits[0].header)
        cutout_low = Cutout2D(
            img_fits[0].data,
            mean_pos,
            sep*img_fac,
            wcs=w,
            copy=True
        )
        new_low_wcs, new_low_shape = find_optimal_celestial_wcs(
            ((cutout_low.data, cutout_low.wcs),), 
            projection='SIN',
            reference=mean_pos
        )


    fig = plt.figure(figsize=(cm(10), cm(24)))

    loc1 = [0.2,0.63,0.8, 0.29]
    loc2 = [0.2, 0.33,0.8, 0.29]
    loc3 = [0.175,0.05,0.78, 0.225]

    img_ax = fig.add_axes(loc1, projection=new_wcs)#cutout.wcs)
    img_low_ax = fig.add_axes(loc2, projection=new_low_wcs)#cutout_low.wcs)
    ax = fig.add_axes(loc3)


    img_ax.imshow(
        cutout.data,
        vmin=-0.001,
        vmax=0.02,
        transform=img_ax.get_transform(cutout.wcs)
    )
    grid_kw = dict(color='white', ls='solid')
    img_ax.coords.grid(True, **grid_kw)
    lon = img_ax.coords[0]
    lon.set_ticklabel_position('t')
    lon.set_axislabel_position('t')
    lon.set_axislabel('Right Ascension')
    img_ax.coords[1].set_axislabel('Declination')
    overlay_box(img_ax, '170-231$\,$MHz', x=0.02, y=0.02)
    if deep_psf is not None:
        add_beam(img_ax, new_wcs, deep_psf, mean_pos, 'upper right', '//')

    fov = 500
    xlim = [-fov, fov]*u.arcsec
    ylim = [-fov, fov]*u.arcsec

    positions = new_wcs.all_world2pix(
            mean_pos.ra + xlim, 
            mean_pos.dec + ylim,
            0
    )
    img_ax.set(
        xlim=positions[0][::-1],
        ylim=positions[1]
    )


    img_low_ax.imshow(
        cutout_low.data,
        vmin=-0.01,
        vmax=0.08,
        transform=img_low_ax.get_transform(cutout_low.wcs)
    )
    img_low_ax.coords.grid(True, **grid_kw)
    lon = img_low_ax.coords[0]
    lon.set_axislabel('Right Ascension')
    img_low_ax.coords[1].set_axislabel('Declination')
    overlay_box(img_low_ax, '72-103$\,$MHz',x=0.02, y=0.02)
    if low_psf is not None:
        add_beam(img_low_ax, new_low_wcs, low_psf, mean_pos, 'upper right', '//')


    fov = 500
    xlim = [-fov, fov]*u.arcsec
    ylim = [-fov, fov]*u.arcsec

    positions = new_low_wcs.all_world2pix(
            mean_pos.ra + xlim, 
            mean_pos.dec + ylim,
            0
    )
    img_low_ax.set(
        xlim=positions[0][::-1],
        ylim=positions[1]
    )

    colours=['red','blue','green']
    markers=['1','x','+']
    tot_flux = []
    tot_fluxerr = []
    for count, (idx, row) in enumerate(isl_df.iterrows()):
        freq, flux, fluxerr = get_freq_flux_err(row, apply_mask=False)
        pos_sky = SkyCoord(row.ref_ra*u.deg, row.ref_dec*u.deg)
        src_str = CoordStr(pos_sky)
        colour = colours[count]
        marker = markers[count]

        tot_flux.append(flux)
        tot_fluxerr.append(fluxerr)

        ax.errorbar(
            freq,
            flux,
            yerr=fluxerr,
            ls='None',
            marker=marker,
            label=src_str.tex,
            color=colour,
            capsize=2
        )

        img_ax.scatter(
            pos_sky.ra, 
            pos_sky.dec, 
            marker=marker,
            s=100,
            transform=img_ax.get_transform('fk5'),
            color=colour, 
            zorder=100
        )

        img_low_ax.scatter(
            pos_sky.ra, 
            pos_sky.dec, 
            marker=marker,
            s=100,
            transform=img_low_ax.get_transform('fk5'),
            color=colour,
            zorder=100
        )

    tot_flux = np.sum(tot_flux, axis=0)
    tot_fluxerr = np.sqrt(np.sum(np.array(tot_fluxerr)**2, axis=0))
    mask = np.isfinite(tot_flux) & np.isfinite(tot_fluxerr)
    ax.errorbar(
        freq[mask],
        tot_flux[mask],
        yerr=tot_fluxerr[mask],
        linestyle='None',
        marker='+',
        label='Total',
        capsize=2
    )

    ax.legend(loc='upper right', prop={'size': 6})
    ax.loglog()
    ax.set(
        xlabel='Frequency (MHz)',
        ylabel='Integrated Flux (Jy)'
    )
    ax.grid(
        which='both',
    )
    ax.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax.yaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))
    
    # fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{idx}.png")
    fig.savefig(f"{OUTDIR}/{idx}.pdf")
    plt.close(fig)    

def wrap(info_dict):
        plot_img_sed(*info_dict)

def find_plot_rats(table_path, image_deep_path, image_low_path, sep=80, deep_psf=None, low_psf=None):
    if not os.path.exists(OUTDIR):
        logger.info(f"Creating {OUTDIR}")
        os.mkdir(OUTDIR)

    sep = sep*u.arcsecond

    logger.info(f"Opening {table_path}")
    df = Table.read(table_path).to_pandas()
    logger.info(f"{len(df)} components read in")

    logger.info('Counting number of unique islands')
    isl_counts = Counter(df.island.values)
    logger.info(f"Total of {len(isl_counts.keys())} counted")

    logger.info('Dropping islands with only one components')
    single_comp = [k for k,v in isl_counts.items() if v == 1]
    mask = df.island.isin(single_comp)
    logger.info(f"Dropping {np.sum(mask)} rows")

    crop_df = df[~mask]
    logger.info(f"Table is {len(crop_df)}")


    
    with Pool(4, maxtasksperchild=24) as pool:
        results =  list(tqdm(pool.imap(wrap, [(isl_id, sub_df, image_deep_path, image_low_path, sep, deep_psf, low_psf) for (isl_id, sub_df) in crop_df.groupby('island')], chunksize=8)))


if __name__ == '__main__':
    parser = ArgumentParser(description='Pull out some ratty SEDs')

    parser.add_argument('table', type=str, help='Path to catalogue table to read')
    parser.add_argument('image_deep', type=str, help='Path to the reference deep source finding image')
    parser.add_argument('image_low', type=str, help='Path to a low frequency image')
    parser.add_argument('--island-sep', type=float, default=80, help='Separation in arcsec between components in an island')
    parser.add_argument('--deep-psf', type=str, default=None, help='Path to the projpsf file for the deep source finding image')
    parser.add_argument('--low-psf', type=str, default=None, help='Path to the projpsf file for the low resolution image')

    args = parser.parse_args()

    find_plot_rats(
        args.table,
        args.image_deep,
        args.image_low,
        sep=args.island_sep,
        deep_psf=args.deep_psf,
        low_psf=args.low_psf,
    )