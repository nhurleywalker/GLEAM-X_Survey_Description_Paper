#!/usr/bin/env bash

python reliability/negatives.py \
	IDR1_subset.fits \
	reliability/XG_170-231MHz_psf_comp_negative_dafix_comp.fits \
	-r ref_ra \
	-d ref_dec \
	-o IDR1_subset_filtered.fits

python catalogue_seds/catalogue_seds.py \
	IDR1_subset_filtered.fits \
	-o IDR1_subset_filtered_SEDs.fits

