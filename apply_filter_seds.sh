#!/usr/bin/env bash

python reliability/negatives.py \
	IDR1_subset.fits \
	reliability/XG_170-231MHz_psf_comp_negative_dafix_comp.fits \
	-r ref_ra \
	-d ref_dec \
	-o test_output_idr1.fits

python catalogue_seds/catalogue_seds.py \
	test_output_idr1.fits \
	-o test_output_idr1_wmodels.fits

