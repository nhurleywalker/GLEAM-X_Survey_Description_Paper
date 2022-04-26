# Catalogue is too big for github

# TODO replace with wget command when on ADC

if [[ ! -e GLEAM_downselect.fits ]]
then
    stilts tpipe /home/tash/Dropbox/Public/GLEAM_EGC_v2.fits \
           cmd='keepcols "RAJ2000 DEJ2000 alpha err_alpha reduced_chi2 int_flux_fit_200 err_int_flux_fit_200 int_flux_wide peak_flux_wide local_rms_wide err_int_flux_wide int_flux_151 err_int_flux_151 local_rms_151"' \
           cmd='select "RAJ2000 >= 4*15 && RAJ2000 <= 13*15 && DEJ2000 >= -32.7 && DEJ2000 <= -20.7 && reduced_chi2 <= 1.93 && int_flux_wide/peak_flux_wide <= 2"' \
           out=GLEAM_downselect.fits
fi
