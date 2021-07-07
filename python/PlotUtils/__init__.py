# This file, and the fact that the other files here are in the subdirectory PlotUtils,
# exist only so that the line 'import PlotUtils' will work in other packages.
#
# See http://docs.python.org/2/tutorial/modules.html#packages if you're curious
# how this works.

# load the C++ objects and bind them into the namespace.
import os
if os.getenv("PLOTUTILSTYPE") == "STANDALONE":
    import LoadPlotUtilsLib
else:
    import LoadPlotUtilsLib
