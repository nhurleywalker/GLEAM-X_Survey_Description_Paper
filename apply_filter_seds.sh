#!/usr/bin/env bash

python reliability/negatives.py \
	IDR1_subset.fits \
	reliability/XG_170-231MHz_psf_comp_negative_dafix_comp.fits \
	-r ref_ra \
	-d ref_dec \
	-o IDR1_subset_filtered.fits

python catalogue_seds/catalogue_seds.py \
	IDR1_subset_filtered.fits \
	-c 2 -o IDR1_subset_filtered_SEDs.fits

python catalogue_names/catalogue_names.py \
    IDR1_subset_filtered_SEDs.fits

python add_ucd/add_ucd.py \
    IDR1_subset_filtered_SEDs_paper.fits \
	--apply