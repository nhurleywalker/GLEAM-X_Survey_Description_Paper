import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 


def pl(nu, norm, alpha):
    return norm * nu **alpha 


def make_ax1(ax1, nu):
    ax1.plot(
        nu,
        pl(nu, 1., -0.8),
    )
    ax1.grid(
        which='both',
    )

    ax1.set(
        xscale='log',
        yscale='log',
        xlabel='Frequency (MHz)',
        ylabel='Flux (Jy',
        title='GLEAM-X$\,$J12345.1$-$9876543.23'
    )

def make_small_ax(ax, nu, xlabel=None, onright=False):
    ax.plot(
        nu,
        pl(nu, 1., -0.8),
    )
    ax.grid(
        which='both',
    )

    ax.set(
        xscale='log',
        yscale='log',
        xlabel=xlabel,
        ylabel='Flux (Jy)',
    )

    if onright:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')



example_nu = np.geomspace(40, 2000, 200)

ax1_loc = (0.1, 0.6, 0.8, 0.35)

ax2_loc = (0.1, 0.325, 0.38, 0.2)
ax3_loc = (0.520, 0.325, 0.38, 0.2)

ax4_loc = (0.1, 0.075, 0.38, 0.2)
ax5_loc = (0.520, 0.075, 0.38, 0.2)


fig = plt.figure(figsize=(10, 8))

ax1 = fig.add_axes(ax1_loc)
make_ax1(ax1, example_nu)

ax2 = fig.add_axes(ax2_loc)
make_small_ax(ax2, example_nu)

ax3 = fig.add_axes(ax3_loc)
make_small_ax(ax3, example_nu, onright=True)

ax4 = fig.add_axes(ax4_loc, sharex=ax2)
make_small_ax(ax4, example_nu, xlabel='Frequency (MHz)')

ax5 = fig.add_axes(ax5_loc, sharex=ax2)
make_small_ax(ax5, example_nu, xlabel='Frequency (MHz)', onright=True)

fig.savefig('test.png')
