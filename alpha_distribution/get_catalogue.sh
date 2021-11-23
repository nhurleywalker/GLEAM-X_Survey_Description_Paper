# Catalogue is too big for github

# TODO replace with wget command when on ADC

if [[ ! -e IDR_v1.1_joined_rescaled_cat_seds.vot ]]
then
    scp nhurleywalker@magnus.pawsey.org.au:/astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/170-231MHz_Reference_Results_Deep_EPS5/IDR_v1.1_joined_rescaled_cat_seds.vot ./
fi

if [[ ! -e IDR_v1.1_joined_rescaled_cat_seds_subset.fits ]]
then
    stilts tpipe IDR_v1.1_joined_rescaled_cat_seds.vot \
           cmd='keepcols "ref_ra ref_dec sp_alpha sp_alpha_err sp_norm sp_norm_err sp_rchi2 int_flux peak_flux local_rms"' \
           cmd='select "ref_ra >= 4*15 && ref_ra <= 13*15 && ref_dec <= -20.7 && ref_dec >= -32.7 && sp_rchi2 <= 1.93 && int_flux/peak_flux <= 2"' \
           out=IDR_v1.1_joined_rescaled_cat_seds_subset.fits
fi
