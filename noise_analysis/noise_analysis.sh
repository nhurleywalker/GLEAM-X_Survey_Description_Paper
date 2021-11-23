#!/bin/bash

# First extract small pieces so we don't run out of memory
ra="10:30:00"
dec="-27:30:00"
imgpix=1250
bkgpix=70

if [[ ! -e XG.fits ]]
then
    getfits -o XG.fits XG_170-231MHz.fits $ra $dec $imgpix $imgpix
fi

exts="bkg rms"
for ext in $exts
do
    if [[ ! -e XG_${ext}.fits ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE XG.fits
    fi
done

# Make a small local catalogue
# TODO only define these numbers once, above, and set up argparse for this python script
if [[ ! -e XG_comp.fits ]]
then
# Make my own catalogue
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --cores 1 \
    --autoload \
    --psf=XG_170-231MHz_projpsf_psf.fits \
    --table=XG.fits \
    XG.fits
fi

# Subtract catalogue from image
/usr/local/bin/AeRes -c XG_comp.fits -f XG.fits -r XG_residual.fits
# Subtract background from image
subtract.py XG.fits XG_bkg.fits XG_bkgsubtracted.fits
# Subtract catalogue from background-subtracted image
/usr/local/bin/AeRes -c XG_comp.fits -f XG_bkgsubtracted.fits -r XG_residual_bkgsubtracted.fits
# Mask instead of subtract
/usr/local/bin/AeRes --mask --sigma=1 -c XG_comp.fits -f XG_bkgsubtracted.fits -r XG_masked_bkgsubtracted.fits
# Create S/N map
div.py XG.fits XG_rms.fits XG_sigma.fits
# Background-subtracted S/N map
div.py XG_bkgsubtracted.fits XG_rms.fits XG_bkgsubtracted_sigma.fits
# Masked, background-subtracted S/N map
div.py XG_masked_bkgsubtracted.fits XG_rms.fits XG_masked_bkgsubtracted_sigma.fits
# Residuals, background-subtracted S/N map
div.py XG_residual_bkgsubtracted.fits XG_rms.fits XG_residual_bkgsubtracted_sigma.fits

