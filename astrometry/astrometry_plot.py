from argparse import ArgumentParser
import logging

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky, match_coordinates_sky
from astropy.table import Table
from scipy.optimize import curve_fit
from tqdm import tqdm
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]}
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def inspect_astrometry(gleam_sky, ref_sky, plot_output=None, info=None):
    info = {} if info is None else info
    
    res = match_coordinates_sky(
                gleam_sky,
                ref_sky,
                nthneighbor=1,
    )
    mask = res[1] < 60*u.arcsecond
    res = [i[mask] for i in res]
    gleam_sub_sky = gleam_sky[mask]
    
    info['reference_matches'] = np.sum(mask)
    info['gleam_sky_count'] = len(gleam_sky)
    
    pos_1 = ref_sky[res[0]]
    frame = pos_1.skyoffset_frame()

    deltas = gleam_sub_sky.transform_to(frame)

    x = deltas.lon.to(u.arcsecond).value
    y = deltas.lat.to(u.arcsecond).value

    info['ra_mean_delta'] = np.mean(x)
    info['dec_mean_delta'] = np.mean(y)
    info['ra_median_delta'] = np.median(x)
    info['dec_median_delta'] = np.median(y)
    info['ra_std_delta'] = np.std(x)
    info['dec_std_delta'] = np.std(y)
    info['mask'] = mask
    
    if plot_output is not None:
        logger.info(f"Creating {plot_output} plot")

        x_lim = [-2,2]
        y_lim = [-2,2]
        
        x_mask = (x_lim[0] < x) & (x <= x_lim[1])
        y_mask = (y_lim[0] < y) & (y <= y_lim[1])
        mask = x_mask & y_mask
        
        x = x[mask]
        y = y[mask]
        
        fig, ax2 = plt.subplots(1,1, figsize=(6,6))

        xy = np.vstack([x, y])

        z = gaussian_kde(xy)(xy)

        idx = np.argsort(z)
        x = x[idx]
        y = y[idx]
        z = z[idx]

        ax2.scatter(
            x,
            y,
            marker='o',
            s=2,
            c=z
        )
        ax2.grid(which='both')
        x_median = np.median(x)
        y_median = np.median(y)

        ax2.axvline(x_median, ls='--', color='black')
        ax2.axhline(y_median, ls='--', color='black')
        
        ax2.set(
            xlim=x_lim,
            ylim=y_lim,
        )

        divider = make_axes_locatable(ax2)
        ax_top = divider.append_axes('top', size='15%', pad='4%', sharex=ax2)
        ax_right = divider.append_axes('right', size='15%', pad='4%', sharey=ax2)

        b = ax_top.hist(
            x, bins=30
        )

        b = ax_right.hist(
            y, bins=30, orientation='horizontal'
        )

        ax_top.axvline(x_median, ls='--', color='black')
        ax_right.axhline(y_median, ls='--', color='black')

        ax_top.xaxis.set_visible(False)
        ax_right.yaxis.set_visible(False)
        
        ax2.set(
            xlabel='RA Offset (arcsecond)',
            ylabel='Dec Offset (arcsecond)',
        )
        
        fig.tight_layout()
        fig.savefig(f"{plot_output}")
        
    return info


def create_astrometry(
    src_cata_path, 
    ref_cata_path, 
    *args,
    src_ra='ra', 
    src_dec='dec', 
    ref_ra='RAJ2000', 
    ref_dec='DEJ2000',
    plot_output=None,
    **kwargs
    ):
    logger.info(f"Loading {src_cata_path}")
    src_cata = Table.read(src_cata_path).to_pandas()

    logger.info(f"Loading reference catalogue {ref_cata_path}")
    ref_cata = Table.read(ref_cata_path).to_pandas()

    logger.info(f'Creating source sky positions using {src_ra} and {src_dec} ')
    src_pos = SkyCoord(
        src_cata[src_ra] * u.deg,
        src_cata[src_dec] * u.deg
    )

    logger.info(f'Creating reference sky positions using {ref_ra} and {ref_dec} ')
    ref_pos = SkyCoord(
        ref_cata[ref_ra] * u.deg,
        ref_cata[ref_dec] * u.deg
    )

    info = inspect_astrometry(src_pos, ref_pos, plot_output=plot_output)

    for k, v in info.items():
        logger.debug(f"{k} {v}")


if __name__ == '__main__':
    parser = ArgumentParser(description='Creates the astrometry plots for the GLEAM-X IDR1 paper')
    parser.add_argument('catalogue', type=str, help='Path to the source catalogue to process')
    parser.add_argument('--ra', type=str, default='ref_ra', help='The name of the RA column in the source catalogue to use')
    parser.add_argument('--dec', type=str, default='ref_dec', help='The name of the DEC column in the source catalogue to use')

    parser.add_argument('ref_catalogue', type=str, help='Path to the source catalogue to process')
    parser.add_argument('--ref-ra', type=str, default='RAJ2000', help='The name of the RA column in the reference catalogue to use')
    parser.add_argument('--ref-dec', type=str, default='DEJ2000', help='The name of the DEC column in the reference catalogue to use')

    parser.add_argument('-p', '--plot-output', default=None, type=str, help='Path to use for the astrometry output plot')

    parser.add_argument('-v','--verbose', default=False, action='store_true', help='Extra output logging')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    create_astrometry(
        args.catalogue,
        args.ref_catalogue,
        src_ra=args.ra,
        src_dec=args.dec,
        ref_ra=args.ref_ra,
        ref_dec=args.ref_dec,
        plot_output=args.plot_output
    )
