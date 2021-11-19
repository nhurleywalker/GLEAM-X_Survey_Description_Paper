# Catalogue is too big for github

# TODO replace with wget command when on ADC
freqs="072-103 103-134 139-170" 

for freq in $freqs
do

    if [[ ! -e IDR1_XG_${freq}MHz_rescaled.fits ]]
    then
        scp nhurleywalker@magnus.pawsey.org.au:/astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/170-231MHz_Reference_Results_Deep_EPS5/IDR1_XG_${freq}MHz_rescaled.fits ./
    fi

done

freq="170-231"
if [[ ! -e IDR1_XG_${freq}MHz_rescaled.fits ]]
then
    scp nhurleywalker@magnus.pawsey.org.au:/astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/XG_${freq}MHz_rescaled.fits ./IDR1_XG_${freq}MHz_rescaled.fits
fi

# Make a Cartesian image that matches our layout
