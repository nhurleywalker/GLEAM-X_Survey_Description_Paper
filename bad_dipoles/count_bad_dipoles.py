#!/usr/bin/env python

import os
from collections import Counter
import requests
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from argparse import ArgumentParser

BASEURL = "http://ws.mwatelescope.org/metadata"


def download(obsid, service="fits"):
    """Downloads a component from the MWA metadata service

    Args:
        obsid (int): Observation ID of interest
        service (str, optional): Metaservice to call upon. Defaults to "fits".

    Returns:
        requests.models.Response: The response from the URL call
    """
    url = f"{BASEURL}/{service}"

    response = requests.get(url, params={"obsid": obsid}, timeout=1.0)

    return response


def parse_service(obsid, service="fits", out_file=None):
    """Call and oarse the response of the metadata service appropriately. 

    Args:
        obsid (int): Observation ID of interest
        service (str, optional): Which aspect of the metadata service to poll. Defaults to "fits".
        out_file (str, optional): Name of the file to write the resonse to. Only relevant for when a metafits file has been pulled from the fits service. 
                                  A default file suffix of 'down.metafits' will be used if this is None when service is fits. Defaults to None.

    Returns:
        [astropy.table.Table,dict]: If service is fits, a table of the tile / input / antenna mappings is returned. Otherwise, a dict of the formated JSON content is returned
    """
    response = download(obsid, service=service)

    if service == "fits":
        out_file = f"{obsid}_down.metafite" if out_file is None else out_file
        with open(out_file, "wb") as outfile:
            for i in response:
                outfile.write(i)

        tab = Table.read(out_file)

        return tab
    else:
        json_res = response.json()

        return json_res


def dead_dipole_antenna_flag(obsid, meta_outfile=None, dead_tolerance=0):
    """Attempt to identify and flag tiles (ANTENNA ins the measurement set) based on the number of dipoles that are flagged. 
    Since the `tileid` is not a one-to-one correspondence to the antenna number, there has to be some mapping. 

    TODO: I think the `con` service call is redundant and information can be pulled from the `obs` service

    Args:
        obsid (int): Observation ID of interest
        meta_outfile (str, optional): Path to save the metafits file to. If None, one is generated.. Defaults to None.
        dead_tolerance (int, optional): Number of dead dipoles before the tile is considered flagged. Defaults to 0.
    """

    tab = parse_service(obsid, service="fits", out_file=meta_outfile)
    obs = parse_service(obsid, service="obs")
    con = parse_service(obsid, service="con")


    bad_dipoles = obs["rfstreams"]["0"]["bad_dipoles"]
    print(bad_dipoles)
    dead_xx_count = Counter([len(bad_dipoles[k][0]) for k in bad_dipoles.keys()])
    dead_yy_count = Counter([len(bad_dipoles[k][1]) for k in bad_dipoles.keys()])
    dead_tot_count = Counter([len(bad_dipoles[k][0])+len(bad_dipoles[k][1]) for k in bad_dipoles.keys()])
    
    print(f'Results for {obsid}')
    print('Dipoles Dead, Occurance, Pol')
    for c, p in zip([dead_xx_count, dead_yy_count, dead_tot_count],['xx','yy','total']):
        for k, i in sorted(c.items(), key=lambda a: a[0]):
            print(f"{k}, {i}, {p}")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Downloads information from the MWA metadata service to identify and flag tiles that contain some number of dead dipoles"
    )
    parser.add_argument("obsid", type=int, help="observation ID to download")
    parser.add_argument(
        "-o", "--output", default=None, type=str, help="Output path to download to"
    )
    parser.add_argument(
        "-d",
        "--dead-tolerance",
        default=0,
        type=int,
        help="Number of acceptable dead dipoles before the tile is flagged. This is summed across the X and Y. ",
    )
    
    args = parser.parse_args()

    dead_dipole_antenna_flag(
        args.obsid,
        meta_outfile=args.output,
        dead_tolerance=args.dead_tolerance,
    )

