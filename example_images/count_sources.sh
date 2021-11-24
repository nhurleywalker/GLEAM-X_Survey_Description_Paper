#!/bin/bash

catalogue=../alpha_distribution/IDR_v1.1_joined_rescaled_cat_seds_subset.fits

#157.5-(3.5/2)
#155.75000000000000000000
#157.5+(3.5/2)
#159.25000000000000000000
#-27.5+(3.5/2)
#-25.75000000000000000000
#-27.5-(3.5/2)
#-29.25000000000000000000

# TODO: synchronise this with the plotting code
ra1=155.75
ra2=159.25
dec1=-25.75
dec2=-29.25


stilts tpipe \
        in=$catalogue \
        cmd="select ((ref_ra<=$ra2)&&(ref_ra>=$ra1)&&(ref_dec<=$dec1)&&(ref_dec>=$dec2)&&(int_flux/local_rms>=5))" \
        cmd='keepcols "ref_ra ref_dec"' \
        out='./number_sources.csv'

n=`wc -l number_sources.csv | awk '{print $1}'`
((n-=1))
echo "$n sources in region"

