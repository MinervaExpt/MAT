#include "MnvApplication.h"

void MAT::Initialize()
{
  //! Enable cintex to get reflex powers with PlotUtils (necessary for I/O and various other ROOT happenings)
#ifndef MNVROOT6
  ROOT::Cintex::Cintex::Enable();
#endif
}
