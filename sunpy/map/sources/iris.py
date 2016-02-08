from __future__ import absolute_import

import numpy as np
import sunpy.cm as cm
from sunpy.map import GenericMap
import matplotlib.colors as colors
from sunpy.visualization import wcsaxes_compat

__all__ = ['SJIMap']

BAD_PIXEL_VALUE = -32768


class SJIMap(GenericMap):
    """
    A 2D IRIS Slit Jaw Imager Map.

    The Interface Region Imaging Spectrograph (IRIS) small explorer spacecraft
    provides simultaneous spectra and images of the photosphere, chromosphere,
    transition region, and corona with 0.33 to 0.4 arcsec spatial resolution,
    2-second temporal resolution and 1 km/s velocity resolution over a
    field-of- view of up to 175 arcsec by 175 arcsec.  IRIS consists of a 19-cm
    UV telescope that feeds a slit-based dual-bandpass imaging spectrograph.

    Slit-jaw images in four different passbands (C ii 1330, Si iv 1400,
    Mg ii k 2796 and Mg ii wing 2830  A) can be taken simultaneously with
    spectral rasters that sample regions up to 130 arcsec by 175 arcsec at a
    variety of spatial samplings (from 0.33 arcsec and up).
    IRIS is sensitive to emission from plasma at temperatures between
    5000 K and 10 MK.

    IRIS was launched into a Sun-synchronous orbit on 27 June 2013.

    References
    ----------
    * `IRIS Mission Page <http://iris.lmsal.com>`_
    * `IRIS Analysis Guide <https://iris.lmsal.com/itn26/itn26.pdf>`_
    * `IRIS Instrument Paper <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0196&file_type=pdf>`_
    * `IRIS FITS Header keywords <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0077&file_type=pdf>`_
    """

    def __init__(self, data, header, **kwargs):
        GenericMap.__init__(self, data, header, **kwargs)
        self.data = np.ma.masked_equal(data, BAD_PIXEL_VALUE)
        if header.get('lvl_num') == 2:
            self.meta['wavelnth'] = header.get('twave1')
            self.meta['detector'] = header.get('instrume')
            self.meta['waveunit'] = "Angstrom"
        if header.get('lvl_num') == 1:
            self.meta['wavelnth'] = int(header.get('img_path').split('_')[1])
            self.meta['waveunit'] = "Angstrom"

        self.meta['detector'] = "SJI"
        self.meta['waveunit'] = "Angstrom"

        self.plot_settings['cmap'] = cm.get_cmap('irissji' + str(int(self.meta['wavelnth'])))
        self.plot_settings['norm'] = colors.PowerNorm(0.4)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an IRIS SJI image"""
        tele = header.get('TELESCOP', '').startswith('IRIS')
        obs = header.get('INSTRUME', '').startswith('SJI')
        level = header.get('lvl_num') == 1
        return tele and obs

    def draw_slit(self, axes=None, **kwargs):
        """Draws the slit location over the SJI observation.

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
        Axes to plot limb on or None to use current axes.

        Returns
        -------
        lines: list
            A list of `matplotlib.axvline` objects that have been plotted.

        Notes
        -----
        keyword arguments are passed onto matplotlib.pyplot.plot
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)


        axes.axvline(self.meta['SLTPX1IX'])
