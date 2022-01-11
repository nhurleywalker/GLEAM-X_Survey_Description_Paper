from argparse import ArgumentParser
import os
from numpy.lib.scimath import power
import logging
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')

from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.table import Table 
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from scipy.stats import chi2


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

REF_NU = 200.
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
    
    return freq, int_flux, err_flux


def power_law(nu, norm, alpha):
    return norm * (nu / REF_NU) ** alpha


def curved_power_law(nu, norm, alpha, q):
    spec_nu = nu / REF_NU
        
    return norm * spec_nu ** alpha * \
            np.exp(q * np.log(spec_nu)**2)

def plot_sed(freq, flux, fluxerr, pl_res, cpl_res, coord_src, outpath):
    fig, ax = plt.subplots(1,1)

    nu = np.geomspace(
        np.min(freq),
        np.max(freq),
        100
    )

    ax.errorbar(
        freq,
        flux,
        yerr=fluxerr,
        ls='None',
        marker='.'
    )

    legend = False

    if pl_res is not None:
        legend = True
        pl_reject =  pl_res['chi2'] > chi2.ppf(0.99, pl_res['dof'])
        pl_params = [pl_res[i] for i in ['norm', 'alpha']]

        ax.plot(
            nu,
            power_law(nu, *pl_params),
            ls='--',
            color='red',
            label=f"Power-law, r$\chi^2$  {pl_res['rchi2']:.3f}, reject {pl_reject}"
        )

    if cpl_res is not None:
        legend = True
        cpl_reject =  cpl_res['chi2'] > chi2.ppf(0.99, cpl_res['dof'])
        cpl_params = [cpl_res[i] for i in ['norm', 'alpha', 'q']]

        ax.plot(
            nu,
            curved_power_law(nu, *cpl_params),
            ls=':',
            color='green',
            label=f'Curved power-law, r$\chi^2$ {cpl_res["rchi2"]:.3f}, reject {cpl_reject}'
        )

    if legend is True:
        ax.legend()

    ax.loglog()
    ax.set(
        xlabel='Frequency (MHz)',
        ylabel='Integrated Flux (Jy)',
        title=coord_src.tex
    )

    fig.tight_layout()
    fig.savefig(f"{outpath}/{coord_src}.png")
    plt.close(fig)


def fit_pl(freq, flux, fluxerr):
    p0 = (np.median(flux), -0.8)

    try:
        fit_res = curve_fit(
            power_law,
            freq,
            flux,
            p0=p0,
            sigma=fluxerr,
            absolute_sigma=True
        )
    except RuntimeError:
        return None

    best_p, covar = fit_res
    err_p = np.sqrt(np.diag(covar))
    dof = len(freq) - 2
    chi2 = np.sum(
        ((flux - power_law(freq, *best_p)) / fluxerr)**2
        )
    rchi2 = chi2 / dof

    return dict(
        norm=best_p[0], 
        alpha=best_p[1],
        norm_err=err_p[0],
        alpha_err = err_p[1],
        chi2=chi2,
        rchi2=rchi2,
        dof=dof
        )


def fit_cpl(freq, flux, fluxerr):
    p0 = (np.median(flux), -0.8, 0)

    try:
        fit_res = curve_fit(
            curved_power_law,
            freq,
            flux,
            p0=p0,
            sigma=fluxerr,
            absolute_sigma=True
        )
    except RuntimeError:
        return None

    best_p, covar = fit_res
    err_p = np.sqrt(np.diag(covar))
    dof = len(freq) - 3
    chi2 = np.sum(
        ((flux - curved_power_law(freq, *best_p)) / fluxerr)**2
        )
    rchi2 = chi2 / dof
    return dict(
        norm=best_p[0], 
        alpha=best_p[1],
        q=best_p[2],
        norm_err=err_p[0],
        alpha_err = err_p[1],
        q_err=err_p[2],
        chi2=chi2,
        rchi2=rchi2,
        dof=dof
        )


def fit_models(info):
    idx = info['idx']
    row = info['row']
    plot = info['plot']
    # idx, row = row

    freq, flux, fluxerr = get_freq_flux_err(row)
    
    if len(freq) < 15:
        return {}

    pl_res_orig = fit_pl(freq, flux, fluxerr)
    cpl_res_orig = fit_cpl(freq, flux, fluxerr)

    pl_res = None if pl_res_orig is None or pl_res_orig['chi2'] > chi2.ppf(0.99, pl_res_orig['dof']) else pl_res_orig
    cpl_res = None if cpl_res_orig is None or cpl_res_orig['chi2'] > chi2.ppf(0.99, cpl_res_orig['dof']) else cpl_res_orig

    if cpl_res is not None and (np.abs(cpl_res['q']) < 0.2 or cpl_res['q']/cpl_res['q_err'] < 3):
        cpl_res = None

    # If both failed, nothing can be done
    if pl_res is None and cpl_res is None:
        return {}
    
    if plot is True:
        coord_src = CoordStr(
            SkyCoord(row.ref_ra*u.deg, row.ref_dec*u.deg)
        )
        plot_sed(freq, flux, fluxerr, pl_res_orig, cpl_res_orig, coord_src, 'SEDs')


    # hate this
    if pl_res is not None and cpl_res is None:
        best_sn, best_res = 'pl', pl_res 
    elif cpl_res is not None and pl_res is None:
        best_sn, best_res = 'cpl', cpl_res 
    elif pl_res['rchi2'] < cpl_res['rchi2']:
        best_sn, best_res = 'pl', pl_res 
    else:
        best_sn, best_res = 'cpl', cpl_res 

    results = {f"{best_sn}_{k}":v for k,v in best_res.items()}
    results['index'] = idx

    return results


def process_catalogue(tab_path, output=None, plot=False, cpus=1, chunksize=16):

    if plot is True:
        if not os.path.exists('SEDs'):
            logger.info('Creating SEDs folder')
            os.mkdir('SEDs')

    df = Table.read(tab_path).to_pandas()
    logger.info(f"Read in {tab_path}, contains {len(df)} rows")

    logger.info(f"Starting fits")
    with Pool(cpus, maxtasksperchild=24) as pool:
        results =  list(tqdm(pool.imap(fit_models, [dict(idx=i, row=r, plot=plot) for i, r in df.iterrows()], chunksize=chunksize)))

    logger.info("Creating results dataframe")
    res_df = pd.DataFrame(results).set_index('index')
    logger.info(f'Result set length is {len(res_df)}')

    logger.info('Combining dataframes')
    comb_df = df.join(res_df)

    pl_count = np.sum(np.isfinite(comb_df.pl_rchi2))
    cpl_count = np.sum(np.isfinite(comb_df.cpl_rchi2))
    no_count = np.sum( ~np.isfinite(comb_df.pl_rchi2) & ~np.isfinite(comb_df.cpl_chi2) )


    logger.info(f"Power law selected: {pl_count}")
    logger.info(f"Curved power law selected: {cpl_count}")
    logger.info(f"No model selected: {no_count}")
    logger.info(f"Total models: {pl_count + cpl_count + no_count}")

    if output is None:
        logger.info('Not writing any saved tables')
    else:
        logger.info('Saving table. converting to astropy table')
        tab = Table.from_pandas(comb_df)

        # Clean up 1e+20 fill values from_pandas() brings in
        cols = [c for c in tab.columns if 'pl_' in c]
        for c in cols:
            logger.info(f"Cleaning column {c}, {np.max(tab[c])}, {np.sum(tab[c].mask)}")
            tab[c][tab[c].mask] = np.nan


        logger.info(f"Writing to {output}")
        tab.write(output, overwrite=True, format='fits')


if __name__ == '__main__':
    parser = ArgumentParser(description='Perform fitting and other todos on catalogue')

    parser.add_argument('table', type=str, help='Path to table to read')
    parser.add_argument('-o','--output', default=None, type=str, help='Output Table to write (using astropy.table.Table)')
    parser.add_argument('-p','--plot', default=False, action='store_true', help='Make plots for SEDs with a successful fit')
    parser.add_argument('-c','--cpus', default=1, type=int, help='Number of CPUs to abuse')
    parser.add_argument('--chunksize', default=16, type=int, help='Chunksize used in the multiprocessing. Bigger numbers can be a lot faster. Be cautious of large numbers when plotting. ')
    args = parser.parse_args()

    process_catalogue(
        args.table,
        output=args.output,
        plot=args.plot,
        cpus=args.cpus
    )
