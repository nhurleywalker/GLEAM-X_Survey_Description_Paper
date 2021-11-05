#!/usr/bin/env python

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

import sys

# TODO: make into argparse
cat = sys.argv[1]
ra = sys.argv[2]
dec = sys.argv[3]
rad = sys.argv[4]
out = sys.argv[5]

hdu = fits.open(cat)
c = SkyCoord(ra, dec, frame="fk5", unit=(u.deg, u.deg))
ccat = SkyCoord(hdu[1].data["ra"], hdu[1].data["dec"], frame="fk5", unit=(u.deg, u.deg))
idx = ccat.separation(c) < float(rad)*u.deg

#idxc, idxcatalog, d2d, d3d = ccat.search_around_sky(c, float(rad)*u.deg)

hdu[1].data = hdu[1].data[idx]

hdu.writeto(out)
