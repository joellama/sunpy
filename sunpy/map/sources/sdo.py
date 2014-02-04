"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from matplotlib import colors
from astropy.units import Quantity

from sunpy.map import GenericMap
from sunpy.cm import cm

__all__ = ['AIAMap', 'HMIMap']

class AIAMap(GenericMap):
    """AIA Image Map definition
    
    References
    ----------
    For a description of AIA headers
    http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf
    """
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        # Fill in some missing info
        self.meta['detector'] = "AIA"
#        self.meta['instrme'] = "AIA"
        
        self._nickname = self.detector
        self._name = self.detector + " " + str(self.measurement)
        
        self.cmap = cm.get_cmap('sdoaia%d' % self.wavelength)
    
    @property
    def observatory(self):
        return self.meta['telescop'].split('/')[0]
        
    @property
    def processing_level(self):
        return self.meta['lvl_num']

    def _get_norm(self):
        """Returns a Normalize object to be used with AIA data"""
        # byte-scaled images have most likely already been scaled
        if self.data.dtype == np.uint8:
            return None

        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        return colors.Normalize(vmin, vmax)
    
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('AIA')
        
class HMIMap(GenericMap):
    """HMI Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        # Mask the array of all nans
        
        self.meta['detector'] = "HMI"
#        self.meta['instrme'] = "HMI"
#        self.meta['obsrvtry'] = "SDO"

        self._name = self.detector + " " + str(self.measurement)
        self._nickname = self.detector
        self.cmap = cm.get_cmap('sdohmimag')
        #the following should really be happening in map base!
        self.data = self.rotate(np.deg2rad(self.rotation_angle.get('y'))).data
        self.data = np.ma.masked_array(self.data,np.isnan(self.data))

        # next set crota2 to zero
        
    def unsigned_magnetic_flux(self):
        """Return the unsigned magnetic flux for the map"""
        # the following value should probably be available in map
        arcsec_to_m = np.sin(np.deg2rad(1/(60.*60.))) * Quantity(self.dsun, 'm')
        pixel_area = arcsec_to_m**2 * self.scale['x'] * self.scale['y']
        xx, yy = self.pixel_to_data()
        angles = np.sqrt(xx**2 + yy**2)/self.rsun_arcseconds * 90.0
        angles = np.ma.masked_array(angles,np.isnan(angles))
        pixel_angles = np.sqrt(xx**2 + yy**2)/self.rsun_arcseconds * 90.0
        factor = 1/np.cos(np.deg2rad(angles))**2
        bflux = pixel_area.to('cm**2') * Quantity(self.data,'gauss') * factor.transpose()
        return np.sum(np.abs(bflux))

    def _get_norm(self):
        """Returns a Normalize object to be used with AIA data"""
        # byte-scaled images have most likely already been scaled
        if self.data.dtype == np.uint8:
            return None
        return colors.Normalize(-300, 300)

    @property
    def measurement(self):
        return self.meta['content'].split(" ")[0].lower()
    
    @property
    def observatory(self):
        return self.meta['telescop'].split('/')[0]    

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume', '').startswith('HMI') 
