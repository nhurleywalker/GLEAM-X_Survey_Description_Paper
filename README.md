# GLEAM-X_Survey_Description_Paper
Plots and data to support the GLEAM-X Survey Description Paper


## Notes on IDR1_subset.fits

The complete priorised catalogue of the region is ~230,000 components, and is ~1.4GB. This `IDR1_subset.fits` in this repository is a catalogue extract created using the following Topcat command:

`ref_ra >= 4*15 && ref_ra <= 13*15 && ref_dec <= -20.7 && ref_dec >= -32.7 && (int_flux/local_rms)>=5`

For completeness, the output `vot` table produced from the `aegean priorised` fitting pipeline may be downloaded from:

<https://cloudstor.aarnet.edu.au/plus/s/h7FtOkZrsaWkwPW>

Be mindful that since it is a `vot` table some codes (particularly `read()` from `astropy.table.Table`) might take a few moments to parse the input.

## apply_filter_seds.sh

A very simple shell script is include to execute the `reliability/negatives.py` and `catalogue_seds/catalogue_seds.py` scripts. 

## git-lfs

The `IDR1_subset.fits` is ~305MB, and has been commited to this repository using `git lfs`. This `git` extension may not be installed by default on some systems. Instructions on the installation and basic use of `git lfs` may be found:

<https://git-lfs.github.com/>