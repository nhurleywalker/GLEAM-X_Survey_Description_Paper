from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import Quadrangle
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

LIMITS = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['14:00:00 -19:48:00','05:00:00 -35:26:00']]


def make_completeness_plots(comp_cube, base_out='Completeness', s_min=-3, s_max=-0.5, s_step=0.1):
    comp_cube_fits = fits.open(comp_cube)

    w = WCS(comp_cube_fits[0].header).celestial

    x, y = np.indices(comp_cube_fits[0].data[0].shape)
    coords = w.wcs_pix2world(y, x, 0)

    mask = (LIMITS[1].ra.deg < coords[0]) & (coords[0] <= LIMITS[0].ra.deg) &\
            (LIMITS[1].dec.deg < coords[1]) & (coords[1] <= LIMITS[0].dec.deg)

    stats = np.nanpercentile(
            comp_cube_fits[0].data[:, mask], 
            [16, 50, 84], 
            axis=1
    )

    print(stats)
    print(np.sum(np.isfinite(comp_cube_fits[0].data[:, mask])))

    s = np.arange(s_min, s_max + 0.0001, s_step)
    sdim = len(s)
    slin = 10 ** (s + 3.0)  # convert to mJy in linear space

    curve = interp1d(
                stats[1],
                slin,
    )

    fig, ax = plt.subplots(1,1, figsize=(4,4))

    ax.errorbar(
        slin,
        stats[1],
        ms=2.0, 
        color="k",
        marker='o',
        linewidth=0.7,
        yerr=[stats[1]-stats[0], stats[2]-stats[1]]
    )

    ax.axhline(
        100,
        ls='--'
    )
    ax.grid(
        'major'
    )

    ax.set(
        xscale='log',
        xlabel='Flux density (mJy)',
        ylabel='Completeness (%)',
    )

    fig.tight_layout()
    fig.savefig(f'{base_out}_curve.pdf')



if __name__ == '__main__':
    parser = ArgumentParser(description='Make some completeness figures up')
    parser.add_argument('completeness_cube', type=str, help='Path to the FITS cube produced by the completeness simulations')
    parser.add_argument('-b','--base-out', type=str, default='Completeness', help='Basename to use when making up output files')

    args = parser.parse_args()

    make_completeness_plots(
        args.completeness_cube,
        base_out=args.base_out
    )