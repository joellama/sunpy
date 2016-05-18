"""
================
AIA Plot Example
================

This is a very simple way to plot a sample AIA image.
"""

from sunpy.data.sample import EIT_195_IMAGE
import sunpy.map
import matplotlib.pyplot as plt
import numpy.ma

###############################################################################
# We now create the Map using the sample data.

eitmap = sunpy.map.Map(EIT_195_IMAGE)
eitmap.mask = numpy.ma.greater(eitmap.data, 100)
cmap = eitmap.plot_settings['cmap']
cmap.set_bad('b', 1.0)
###############################################################################
# Now we do a quick plot.

eitmap.peek()
plt.show()
