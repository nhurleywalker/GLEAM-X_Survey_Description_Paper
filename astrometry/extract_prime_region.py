#!/usr/bin/env python

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

import sys

# TODO: make into argparse
cat = sys.argv[1]
out = sys.argv[2]

# Nights in IDR1   RA start    RA end
# 2018-02-04         50           220
# 2018-02-09         40           215
# 2018-02-20         60           240
# 2018-03-05         85           270
# So at maximum 50 to 240
# Also need to subtract ~10 deg from each end because of PB cutoff
# 60 to 215
# Visually, image quality is much lower past 190, probably because of Cen A
# So let's go with 60 to 195 = 4h to 13h

ramin = 60
ramax = 195
decmin = -36.7
decmax = -16.7

hdu = fits.open(cat)
idx = np.where(hdu[1].data["ra"] > ramin)
idx = np.intersect1d(idx, np.where(hdu[1].data["ra"] < ramax))
idx = np.intersect1d(idx, np.where(hdu[1].data["dec"] > decmin))
idx = np.intersect1d(idx, np.where(hdu[1].data["dec"] < decmax))

hdu[1].data = hdu[1].data[idx]

hdu.writeto(out)
