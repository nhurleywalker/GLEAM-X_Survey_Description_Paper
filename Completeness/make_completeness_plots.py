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
from mpl_toolkits.axes_grid1.mpl_axes import Axes

LIMITS = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['14:00:00 -19:48:00','05:00:00 -35:26:00']]


def make_curve_plot(stats, flux_levels, base_out):
    curve = interp1d(
                stats[1],
                flux_levels,
    )

    fig, ax = plt.subplots(1,1, figsize=(4,4))

    ax.errorbar(
        flux_levels,
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


def overlay_box(ax, text, x=0.02, y=0.125):
    ax.text(
        x, 
        y,
        text, 
        transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.4, edgecolor='black', boxstyle='round,pad=0.25')
    )


def make_spatial_plot(comp_cube, flux_levels, w, base_out, cmap='inferno'):
    
    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(4,1,1, projection=w)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('top', pad=0, size='5%', axes_class=Axes)

    ax1.imshow(
        comp_cube[0].data[6],
        vmin=0,
        vmax=100, cmap='inferno'
    )
    ax1.set(
        xlim=[149, 284],
        ylim=[50,70],
        ylabel='Dec'
    )
    overlay_box(ax1, f"{flux_levels[6]:.2f} mJy")
    overlay_box(ax1, "(a)", y=0.75)

    lon = ax1.coords[0]
    lon.set_ticklabel_visible(False)
    lon.set_axislabel('')


    ax2 = fig.add_subplot(4,1,2, projection=w, sharex=ax1)
    # ax2 = divider.append_axes('bottom', pad='3%', size='100%')
    ax2.imshow(
        comp_cube[0].data[9],
        vmin=0,
        vmax=100, 
        cmap='inferno'
    )
    ax2.set(
        xlim=[149, 284],
        ylim=[50,70],
        ylabel='Dec'
    )
    lon = ax2.coords[0]
    lon.set_ticklabel_visible(False)
    lon.set_axislabel('')
    overlay_box(ax2, f"{flux_levels[9]:.2f} mJy")
    overlay_box(ax2, "(b)", y=0.75)
    

    ax3 = fig.add_subplot(4,1,3, projection=w)
    # ax3 = divider.append_axes('bottom', pad='3%', size='100%')
    ax3.imshow(
        comp_cube[0].data[12],
        vmin=0,
        vmax=100, 
        cmap='inferno'
    )
    ax3.set(
        xlim=[149, 284],
        ylim=[50,70],
        ylabel='Dec'
    )
    lon = ax3.coords[0]
    lon.set_ticklabel_visible(False)
    lon.set_axislabel('')
    overlay_box(ax3, f"{flux_levels[12]:.2f} mJy")
    overlay_box(ax3, "(c)", y=0.75)


    ax4 = fig.add_subplot(4,1,4, projection=w)
    # divider = make_axes_locatable(ax3)
    # ax4 = divider.append_axes('bottom', pad='3%', size='100%')
    cim = ax4.imshow(
        comp_cube[0].data[15],
        vmin=0,
        vmax=100, 
        cmap='inferno'
    )
    ax4.set(
        xlim=[149, 284],
        ylim=[50,70],
        xlabel='RA',
        ylabel='Dec'
    )
    overlay_box(ax4, f"{flux_levels[15]:.2f} mJy")
    overlay_box(ax4, "(d)", y=0.75)


    cbar = fig.colorbar(cim, cax=cax, label='Completeness (%)', orientation='horizontal')

    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')


    fig.savefig(f"{base_out}_spatial.png")
    fig.savefig(f"{base_out}_spatial.pdf")


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

    s = np.arange(s_min, s_max + 0.0001, s_step)
    sdim = len(s)
    slin = 10 ** (s + 3.0)  # convert to mJy in linear space

    print(stats)
    print(np.sum(np.isfinite(comp_cube_fits[0].data[:, mask])))


    make_curve_plot(stats, slin, base_out)
    make_spatial_plot(comp_cube_fits, slin, w, base_out)


if __name__ == '__main__':
    parser = ArgumentParser(description='Make some completeness figures up')
    parser.add_argument('completeness_cube', type=str, help='Path to the FITS cube produced by the completeness simulations')
    parser.add_argument('-b','--base-out', type=str, default='Completeness', help='Basename to use when making up output files')

    args = parser.parse_args()

    make_completeness_plots(
        args.completeness_cube,
        base_out=args.base_out
    )