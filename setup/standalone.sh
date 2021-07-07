# needed for PlotUtils to build and work standalone if made with OSXMakefile

export LD_LIBRARY_PATH=${PLOTUTILSROOT}:${LD_LIBRARY_PATH}
export PYTHONPATH=${PLOTUTILSROOT}/python:${PYTHONPATH}
export PLOTUTILSTYPE="STANDALONE"
export PLOTUTILSVERSION="ROOT6"