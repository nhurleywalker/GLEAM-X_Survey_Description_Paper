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


def make_plot(comp_path, coeffs_path, out_path=None):
    df = Table.read(comp_path).to_pandas()
    coeffs = np.loadtxt(coeffs_path)
    poly = np.poly1d(coeffs)

    dec_min = np.min(df.Dec)
    dec_max = np.max(df.Dec)
    
    dec_range = np.linspace(
        dec_min,
        dec_max,
        100
    )

    fig, ax1 = plt.subplots(1,1,figsize=(4,4))

    divider = make_axes_locatable(ax1)
    axt = divider.append_axes('top', size='100%', pad=0)
    axr = divider.append_axes('right', size='25%', pad=0, sharey=ax1)

    density_scatter(
        df.Dec,
        df.log10ratio,
        axt,
        alpha=0.003,
        marker='.',
        edgecolor='None'
    )
    axt.plot(
        dec_range,
        poly(dec_range),
        color='black',
        ls='--',
        lw=2
    )
    axt.xaxis.set_ticklabels([])

    density_scatter(
        df.Dec,
        df.log10ratio_after_full_cor,
        ax1,
        alpha=0.003,
        marker='.',
        edgecolor='None'
    )

    axr.hist(
        df.log10ratio_after_full_cor,
        bins=20,
        orientation='horizontal',
        histtype='step'
    )
    axr.yaxis.set_tick_params(labelleft=False)

    ax1.set(
        xlabel='Declination (Deg)',
        xlim=[dec_min, dec_max],
        ylabel='Corrected Ratio'
    )
    ax1.grid(axis='both')
    axt.set(
        xlim=[dec_min, dec_max],
        ylabel='Log Ratio'
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

    args = parser.parse_args()

    make_plot(
        args.comp_path,
        args.coeffs,
        out_path=args.out_path
    )