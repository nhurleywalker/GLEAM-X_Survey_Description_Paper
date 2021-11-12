from argparse import ArgumentParser

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
from astropy.table import Table 
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 10}
)

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram

    https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
    """
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )

    
    return ax

def weighted_plot(x, y, weight, ax, *args, **kwargs):

    order = np.argsort(weight)
    x = x[order]
    y = y[order]
    w = weight[order]

    SNR = np.log10(w)

    ax.scatter(
        x,
        y,
        marker='.',
        c=SNR,
        cmap='Greys'
    )

    return ax


def make_plot(comp_path, coeffs_path, out_path=None, nsrcs=None):
    df = Table.read(comp_path).to_pandas()
    
    if nsrcs is not None:
        df['weight'] = df['flux'] / df['local_rms']
        df = df.sort_values('weight', ascending=False)
        df = df[:nsrcs]

    coeffs = np.loadtxt(coeffs_path)
    poly_order = len(coeffs) - 1

    dec_min = np.min(df.Dec)
    dec_max = np.max(df.Dec)

    ratio_min, ratio_max = -0.1, 0.1
    
    dec_range = np.linspace(
        dec_min,
        dec_max,
        100
    )

    dec_corr = np.zeros_like(dec_range)
    for i in range(0, poly_order + 1):
        dec_corr += coeffs[i] * pow(dec_range, poly_order - i)

    fig, ax1 = plt.subplots(1,1,figsize=(4,4))

    divider = make_axes_locatable(ax1)
    axt = divider.append_axes('top', size='100%', pad=0)
    axr = divider.append_axes('right', size='25%', pad=0, sharey=ax1)

    weighted_plot(
        df.Dec.values,
        df.log10ratio.values,
        (df.flux / df.local_rms).values,
        axt,
        marker='.',
        edgecolor='None'
    )
    axt.plot(
        dec_range,
        dec_corr,
        ls='-',
        lw=2
    )
    axt.xaxis.set_ticklabels([])

    weighted_plot(
        df.Dec.values,
        df.log10ratio_after_full_cor.values,
        (df.flux / df.local_rms).values,
        ax1,
        marker='.',
        edgecolor='None'
    )

    axr.hist(
        df.log10ratio_after_full_cor,
        bins=25,
        range=[ratio_min, ratio_max],
        orientation='horizontal',
        histtype='step'
    )
    axr.yaxis.set_tick_params(labelleft=False)

    ax1.set(
        xlabel='Declination ($^\circ$)',
        xlim=[dec_min, dec_max],
        ylim=[ratio_min, ratio_max],
        ylabel='$\log_{10} R$ (corrected)'
    )
    ax1.grid(axis='both')
    axt.set(
        xlim=[dec_min, dec_max],
        ylim=[ratio_min, ratio_max],
        ylabel='$\log_{10} R$ (uncorrected)'
    )
    axt.grid(axis='both')

    axr.set(
        xlabel='Counts'
    )

    fig.tight_layout()

    if out_path is None:
        plt.show()
    else:
        fig.savefig(out_path)


if __name__ == '__main__':
    parser = ArgumentParser(description='Make a simple plot for the rescaling magic')
    parser.add_argument('comp_path', type=str, help='Path to rescale component catalogue output')
    parser.add_argument('coeffs', type=str, help='Path to numpy polunomial coefficents')
    parser.add_argument('-o', '--out-path', default=None, type=str, help='Output file to save plot as')
    parser.add_argument('-n', '--nsrcs', type=int, default=10000, help='Numberof sources to plot (nrightes first)')

    args = parser.parse_args()

    make_plot(
        args.comp_path,
        args.coeffs,
        out_path=args.out_path,
        nsrcs=args.nsrcs
    )
