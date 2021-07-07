#ifndef MNVCOLORS_H
#define MNVCOLORS_H

#ifndef __CINT__
#include <cassert>
#include <iostream>
#include <vector>

#include "TClass.h"
#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TList.h"
#include "TPad.h"
#include "TString.h"

namespace MnvColors {
enum EColorPalettes {
  kAlphabetPalette,
  kKellyPalette,
  k36Palette,
  kGlasbeyPalette,
  kBrewerSet1Palette,
  kBrewerDark2Palette,
  kTolBrightPalette,
  kTolMutedPalette,
  kTolLightPalette,
  kOkabeItoPalette,
  kOkabeItoLightPalette,
  kOkabeItoDarkPalette,
  kHeliumPalette,
  kNPalettes
};
const std::vector<int>& GetColors(int palette = kOkabeItoPalette);
void AutoColorHists(TPad* pad, int palette = kOkabeItoPalette);
std::vector<TH1*> GetPadHists(TPad* pad);
const char* GetPaletteName(EColorPalettes e = kNPalettes);
}  // namespace MnvColors

#endif  // __CINT__
#endif  // MNVCOLORS_H
