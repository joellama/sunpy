from __future__ import absolute_import

import numpy as np
import sunpy.cm as cm
from sunpy.map import GenericMap

__all__ = ['IRISMap']

BAD_PIXEL_VALUE = -32768

class IRISMap(GenericMap):
    """
    A 2D IRIS Map
    """

    def __init__(self, data, header, **kwargs):
        GenericMap.__init__(self, data, header, **kwargs)
        self.data = np.ma.masked_equal(self.data, BAD_PIXEL_VALUE)
        self.meta['waveunit'] = 'Angstrom'
        self.meta['wavelnth'] = int(self.meta.get('img_path').split('_')[1])
        self.cmap = cm.get_cmap('irissji' + str(self.meta['wavelnth']))

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        tele = header.get('TELESCOP', '').startswith('IRIS')
        obs = header.get('INSTRUME', '').startswith('SJI')
        return tele and obs
