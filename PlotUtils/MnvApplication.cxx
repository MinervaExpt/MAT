#include "MnvApplication.h"

void PlotUtils::Initialize()
{
  //! Enable cintex to get reflex powers with PlotUtils (necessary for I/O and various other ROOT happenings)
#ifndef MNVROOT6
  ROOT::Cintex::Cintex::Enable();
#endif
}
