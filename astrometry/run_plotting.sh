export PYTHONPATH=
#python3 astrometry_plot.py ../noise_analysis/XG_170-231MHz_psf_comp.fits ~/Programs/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ra --dec dec --plot-output astrometry_for_talks.png --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux
python3 extract_prime_region.py ../noise_analysis/XG_170-231MHz_psf_comp.fits XG_170-231MHz_prime_comp.fits
python3 astrometry_plot.py XG_170-231MHz_prime_comp.fits ~/Programs/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ra --dec dec --plot-output astrometry.pdf --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux --latex
