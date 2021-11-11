#list=`ls XG_D-27_20180220_cenchan_??.txt XG_D-27_20180220_cenchan_???.txt` ; for file in $list ; do obs=`tail -1 $file` ; scp $obs/${obs}_deep-MFS-image-pb_warp_comp.fits tash@glados.cira.curtin.edu.au:/media/data/MWA/GLEAMX/GLEAM-X_Survey_Description_Paper/astrometry/

for file in *_deep-MFS-image-pb_warp_comp.fits
do
    obs=${file:0:10}
#    stilts tmatch2 in1=${file} in2=/home/tash/Programs/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits \
#           values1="ra dec" values2="RAJ2000 DEJ2000" \
#           matcher=sky params=15 out=${obs}_after_warp_xm.fits
    stilts tpipe in=${obs}_after_warp_xm.fits cmd="addcol del_ra 3600*(ra-RAJ2000)" cmd="addcol del_dec 3600*(dec-DEJ2000)" \
           cmd="keepcols 'del_ra del_dec'" \
           omode=stats
done
