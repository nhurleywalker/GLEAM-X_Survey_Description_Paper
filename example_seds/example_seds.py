import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

def pl(nu, norm, alpha):
    return norm * nu **alpha 


def make_ax1(ax1, nu):
    ax1.plot(
        nu,
        pl(nu, 100., -0.8),
    )
    ax1.grid(
        which='both',
    )

    ax1.set(
        xscale='log',
        yscale='log',
        xlabel='Frequency (MHz)',
        ylabel='Flux density (Jy)',
        title='GLEAM-X$\,$J12345.1$-$9876543.23'
    )
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_minor_formatter(FormatStrFormatter('%3.1f'))
#    ax1.xaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))

def make_small_ax(ax, nu, xlabel=None, onright=False):
    ax.plot(
        nu,
        pl(nu, 100., -0.8),
    )
    ax.grid(
        which='both',
    )

    ax.set(
        xscale='log',
        yscale='log',
        ylabel='Flux density (Jy)',
    )
    
    if onright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')

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


example_nu_large = np.geomspace(65, 2000, 200)
example_nu = np.geomspace(65, 270, 200)

ax1_loc = (0.1, 0.6, 0.8, 0.35)
ax2_loc = (0.1, 0.315, 0.395, 0.2)
ax3_loc = (0.505, 0.315, 0.395, 0.2)
ax4_loc = (0.1, 0.075, 0.395, 0.23)
ax5_loc = (0.505, 0.075, 0.395, 0.23)


fig = plt.figure(figsize=(7, 5))

ax1 = fig.add_axes(ax1_loc)
make_ax1(ax1, example_nu_large)

ax2 = fig.add_axes(ax2_loc)
ax3 = fig.add_axes(ax3_loc)
ax4 = fig.add_axes(ax4_loc)#, sharex=ax2)
ax5 = fig.add_axes(ax5_loc)#, sharex=ax3)

print('Top Left')
make_small_ax(ax2, example_nu)

print('Top Right')
make_small_ax(ax3, example_nu, onright=True)

print('Bottom Left')
make_small_ax(ax4, example_nu, xlabel='Frequency (MHz)')

print('Bottom Right')
make_small_ax(ax5, example_nu, xlabel='Frequency (MHz)', onright=True)

fig.savefig('test.png')
