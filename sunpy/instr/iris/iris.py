"""
Some very beta tools for IRIS
"""

import sunpy.io
import sunpy.time
import sunpy.map
from astropy.utils.data import download_file
import os.path
from astropy.io import fits
from sunpy.util.net import url_exists
from shutil import move
from os.path import basename
from sunpy.time import parse_time
from datetime import timedelta

__all__ = ['SJI_to_cube']

def SJI_to_cube(filename):
    """
    Read a SJI file and return a MapCube

    .. warning::
        This function is a very early beta and is not stable. Further work is
        on going to improve SunPy IRIS support.

    Parameters
    ----------
    filename: string
        File to read

    Returns
    -------
    iris_cube: sunpy.map.MapCube
        A map cube of the SJI sequence
    """

    hdus = sunpy.io.read_file(filename)
    #Get the time delta
    time_range = sunpy.time.TimeRange(hdus[hdu][1]['STARTOBS'], hdus[hdu][1]['ENDOBS'])
    splits = time_range.split(hdus[hdu][0].shape[0])

    if not stop:
        stop = len(splits)

    headers = [hdus[hdu][1]]*(stop-start)
    datas = hdus[hdu][0][start:stop]

    #Make the cube:
    iris_cube = sunpy.map.Map(zip(datas,headers),cube=True)
    #Set the date/time
    for i,m in enumerate(iris_cube):
        m.meta['DATE-OBS'] = splits[i].center.isoformat()

    return iris_cube

def get_level1_file(level2_file):
    """Given a IRIS SJI Level 2 file download the associated level 1 individual
    SJI files"""
    base_url = 'http://www.lmsal.com/solarsoft/'
    f = fits.open(level2_file)

    filenames = [a[3][1:] for a in f[2].data]

    for full_file_name in filenames:
#        if url_exists(os.path.join(base_url, full_file_name)):
        l = download_file(os.path.join(base_url, full_file_name))
        real_name = basename(full_file_name)
        print(real_name)
        move(l, os.path.join('/Users/schriste/Desktop/', real_name))

def make_times(level2_file):
    """Return the times of each of the exposures in a level 2 SJI file

    ..note:

    The results were checked with the original level 1 files and agreed to within
    0.1 seconds"""
    iris = fits.open(level2_file)
    times = [parse_time(iris[0].header['STARTOBS']) + timedelta(seconds=dt) for dt in iris[1].data[:,0]]
    return times
