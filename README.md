MINERvA Analysis Toolkit - An error analysis toolkit for particle physics.


[Arxiv 2103.08677][1]


Directory Structure:
  - PlotUtils - Macros and classes to help with ROOT analysis and plotting.


Setup on FNAL machines:
  source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
  setup root v6_10_04d -q e14:prof
  setup cmake v3_7_1
  git clone https://github.com/MinervaExpt/MAT.git
  mkdir -p opt/MATBuild
  cd opt/MATBuild
  cmake ../../MAT -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release
  make install
  cd ../../
  source opt/bin/setup_MAT.sh


[1]: https://arxiv.org/abs/2103.08677
