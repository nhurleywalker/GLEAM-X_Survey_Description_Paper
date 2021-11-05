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

if [[ ! -d XG.im ]]
then
    fits op=xyin in=XG.fits out=XG.im
fi

exts="bkg rms"
for ext in $exts
do
    if [[ ! -e XG_${ext}.fits ]]
    then
    # Get a slightly oversize piece
        getfits -o XG_temp_$ext.fits XG_170-231MHz_$ext.fits $ra $dec $bkgpix $bkgpix
    # Can't use SR6 because it doesn't understand that we've cropped the piece out, so use miriad
        fits op=xyin in=XG_temp_$ext.fits out=XG_temp_$ext.im
        regrid tin=XG.im in=XG_temp_$ext.im out=XG_temp_${ext}_rg.im
        fits op=xyout in=XG_temp_${ext}_rg.im out=XG_${ext}.fits
    # Delete intermediate files
        rm -rf XG_temp_$ext.fits XG_temp_$ext.im XG_temp_${ext}_rg.im
    fi
done

# Make a small local catalogue
# TODO only define these numbers once, above, and set up argparse for this python script
if [[ ! -e XG_subset.fits ]]
then
    python3 extract_local_subset.py XG_170-231MHz_psf_comp.fits 157.5 -27.5 10 XG_subset.fits
fi

# Subtract catalogue from image
/usr/local/bin/AeRes -c XG_subset.fits -f XG.fits -r XG_residual.fits
# Subtract background from image
subtract.py XG.fits XG_bkg.fits XG_bkgsubtracted.fits
# Subtract catalogue from background-subtracted image
/usr/local/bin/AeRes -c XG_subset.fits -f XG_bkgsubtracted.fits -r XG_residual_bkgsubtracted.fits
# Mask instead of subtract
/usr/local/bin/AeRes --mask --sigma=1 -c XG_subset.fits -f XG_bkgsubtracted.fits -r XG_masked_bkgsubtracted.fits
# Create S/N map
div.py XG.fits XG_rms.fits XG_sigma.fits
# Background-subtracted S/N map
div.py XG_bkgsubtracted.fits XG_rms.fits XG_bkgsubtracted_sigma.fits
# Masked, background-subtracted S/N map
div.py XG_masked_bkgsubtracted.fits XG_rms.fits XG_masked_bkgsubtracted_sigma.fits
# Residuals, background-subtracted S/N map
div.py XG_residual_bkgsubtracted.fits XG_rms.fits XG_residual_bkgsubtracted_sigma.fits

#getfits $filedir/Week2_white_lownoise_ddmod_rescaled.fits -o image.fits 9151-12142 8135-11676
#getfits $filedir/Week2_white_lownoise_ddmod_rescaled_blank.fits -o blank_subtracted.fits 9151-12142 8135-11676
#getfits $filedir/xp_bkg.fits -o bkg.fits 9151-12142 8135-11676
#getfits $filedir/xp_rms.fits -o rms.fits 9151-12142 8135-11676
#getfits $filedir/subtracted.fits -o subtracted.fits 9151-12142 8135-11676
# Use hacked version
#AeRes.py -c $filedir/Week2_white_lownoise_ddmod_rescaled_comp.vot -f image.fits -r blanked_image.fits
#AeRes.py -c $filedir/Week2_white_lownoise_ddmod_rescaled_comp.vot -f subtracted.fits -r blanked_subtracted.fits
#AeRes.py -c $filedir/Week2_white_lownoise_ddmod_rescaled_comp.vot -f subtracted.fits -r blanked_subtracted_fc2.fits
# Use unhacked version
#AeRes.py -c $filedir/Week2_white_lownoise_ddmod_rescaled_comp.vot -f subtracted.fits -r residual_subtracted.fits

# Crop the edges to remove badly-subtracted sources
#getfits -o temp.fits 10-2982 10-3532 blanked_subtracted.fits
#mv temp.fits blanked_subtracted.fits
#getfits -o temp.fits 10-2982 10-3532 blanked_subtracted_fc2.fits
#mv temp.fits blanked_subtracted_fc2.fits
#getfits -o temp.fits 10-2982 10-3532 blanked_image.fits
#mv temp.fits blanked_image.fits
#getfits -o temp.fits 10-2982 10-3532 blank_subtracted.fits
#mv temp.fits blank_subtracted.fits
#getfits -o temp.fits 10-2982 10-3532 subtracted.fits
#mv temp.fits subtracted.fits
#getfits -o temp.fits 10-2982 10-3532 residual_subtracted.fits
#mv temp.fits residual_subtracted.fits
#getfits -o temp.fits 10-2982 10-3532 image.fits
#mv temp.fits image.fits
#getfits -o temp.fits 10-2982 10-3532 rms.fits
#mv temp.fits rms.fits
#getfits -o temp.fits 10-2982 10-3532 bkg.fits
#mv temp.fits bkg.fits

