#ifndef MNVCOLORS_CXX
#define MNVCOLORS_CXX

#include "MnvColors.h"

//======================================================================
const std::vector<int>& MnvColors::GetColors(int palette) {
  const int AlphabetColors[] = {
      TColor::GetColor(25, 25, 25),     // ebony
      TColor::GetColor(240, 163, 255),  // amethyst
      TColor::GetColor(255, 0, 16),     // red
      TColor::GetColor(128, 128, 128),  // iron
      TColor::GetColor(25, 164, 5),     // orpiment
      TColor::GetColor(0, 153, 143),    // turquoise
      TColor::GetColor(0, 51, 128),     // navy
      TColor::GetColor(94, 241, 242),   // sky
      TColor::GetColor(157, 204, 0),    // lime
      TColor::GetColor(153, 63, 0),     // caramel
      TColor::GetColor(224, 255, 102),  // uranium
      TColor::GetColor(76, 0, 92),      // damson
      TColor::GetColor(255, 168, 187),  // pink
      TColor::GetColor(194, 0, 136),    // mallow
      TColor::GetColor(255, 204, 153),  // honeydew
      TColor::GetColor(43, 206, 72),    // green
      TColor::GetColor(0, 117, 220),    // blue
      TColor::GetColor(0, 92, 49),      // forest
      TColor::GetColor(255, 225, 0),    // yellow
      TColor::GetColor(153, 0, 0),      // wine
      TColor::GetColor(143, 124, 0),    // khaki
      TColor::GetColor(16, 10, 255),    // violet
      TColor::GetColor(255, 255, 128),  // xanthin
      TColor::GetColor(255, 80, 0),     // zinnia
      TColor::GetColor(148, 255, 181),  // jade
      TColor::GetColor(66, 102, 0)      // quagmire
  };

  const int KellyColors[] = {
      // TColor::GetColor("#f2f3f4"), // white,
      TColor::GetColor("#222222"),  // black,
      TColor::GetColor("#f3c300"),  // yellow,
      TColor::GetColor("#875692"),  // purple,
      TColor::GetColor("#f38400"),  // orange,
      TColor::GetColor("#a1caf1"),  // lightblue,
      TColor::GetColor("#be0032"),  // red,
      TColor::GetColor("#c2b280"),  // buff,
      TColor::GetColor("#848482"),  // gray,
      TColor::GetColor("#008856"),  // green,
      TColor::GetColor("#e68fac"),  // purplishpink,
      TColor::GetColor("#0067a5"),  // blue,
      TColor::GetColor("#f99379"),  // yellowishpink,
      TColor::GetColor("#604e97"),  // violet,
      TColor::GetColor("#f6a600"),  // orangeyellow,
      TColor::GetColor("#b3446c"),  // purplishred,
      TColor::GetColor("#dcd300"),  // greenishyellow,
      TColor::GetColor("#882d17"),  // reddishbrown,
      TColor::GetColor("#8db600"),  // yellowgreen,
      TColor::GetColor("#654522"),  // yellowishbrown,
      TColor::GetColor("#e25822"),  // reddishorange,
      TColor::GetColor("#2b3d26")   // olivegreen
  };

  const int Palette36Colors[] = {
      TColor::GetColor("#5A5156"),
      // TColor::GetColor("#E4E1E3"), // much too light
      TColor::GetColor("#F6222E"),
      TColor::GetColor("#FE00FA"),
      TColor::GetColor("#16FF32"),
      TColor::GetColor("#3283FE"),
      TColor::GetColor("#FEAF16"),
      TColor::GetColor("#B00068"),
      TColor::GetColor("#1CFFCE"),
      TColor::GetColor("#90AD1C"),
      TColor::GetColor("#2ED9FF"),
      TColor::GetColor("#DEA0FD"),
      TColor::GetColor("#AA0DFE"),
      TColor::GetColor("#F8A19F"),
      TColor::GetColor("#325A9B"),
      TColor::GetColor("#C4451C"),
      TColor::GetColor("#1C8356"),
      TColor::GetColor("#85660D"),
      TColor::GetColor("#B10DA1"),
      TColor::GetColor("#FBE426"),
      TColor::GetColor("#1CBE4F"),
      TColor::GetColor("#FA0087"),
      TColor::GetColor("#FC1CBF"),
      TColor::GetColor("#F7E1A0"),
      TColor::GetColor("#C075A6"),
      TColor::GetColor("#782AB6"),
      TColor::GetColor("#AAF400"),
      TColor::GetColor("#BDCDFF"),
      TColor::GetColor("#822E1C"),
      TColor::GetColor("#B5EFB5"),
      TColor::GetColor("#7ED7D1"),
      TColor::GetColor("#1C7F93"),
      TColor::GetColor("#D85FF7"),
      TColor::GetColor("#683B79"),
      TColor::GetColor("#66B0FF"),
      TColor::GetColor("#3B00FB"),
  };

  const int GlasbeyColors[] = {
      TColor::GetColor("#0000FF"), TColor::GetColor("#FF0000"),
      TColor::GetColor("#00FF00"), TColor::GetColor("#000033"),
      TColor::GetColor("#FF00B6"), TColor::GetColor("#005300"),
      TColor::GetColor("#FFD300"), TColor::GetColor("#009FFF"),
      TColor::GetColor("#9A4D42"), TColor::GetColor("#00FFBE"),
      TColor::GetColor("#783FC1"), TColor::GetColor("#1F9698"),
      TColor::GetColor("#FFACFD"), TColor::GetColor("#B1CC71"),
      TColor::GetColor("#F1085C"), TColor::GetColor("#FE8F42"),
      TColor::GetColor("#DD00FF"), TColor::GetColor("#201A01"),
      TColor::GetColor("#720055"), TColor::GetColor("#766C95"),
      TColor::GetColor("#02AD24"), TColor::GetColor("#C8FF00"),
      TColor::GetColor("#886C00"), TColor::GetColor("#FFB79F"),
      TColor::GetColor("#858567"), TColor::GetColor("#A10300"),
      TColor::GetColor("#14F9FF"), TColor::GetColor("#00479E"),
      TColor::GetColor("#DC5E93"), TColor::GetColor("#93D4FF"),
      TColor::GetColor("#004CFF")};

  const int BrewerSet1Colors[] = {
      TColor::GetColor("#e41a1c"), TColor::GetColor("#377eb8"),
      TColor::GetColor("#4daf4a"), TColor::GetColor("#984ea3"),
      TColor::GetColor("#ff7f00"), TColor::GetColor("#ffff33"),
      TColor::GetColor("#a65628"), TColor::GetColor("#f781bf")};

  const int BrewerDark2Colors[] = {
      TColor::GetColor(27, 158, 119),  TColor::GetColor(217, 95, 2),
      TColor::GetColor(117, 112, 179), TColor::GetColor(231, 41, 138),
      TColor::GetColor(102, 166, 30),  TColor::GetColor(230, 171, 2),
      TColor::GetColor(166, 118, 29),  TColor::GetColor(102, 102, 102)};

  // Colorblind-friendly palettes
  // DOI: 10.5281/zenodo.3381072
  const int TolBrightColors[] = {
      TColor::GetColor("#EE6677"), TColor::GetColor("#228833"),
      TColor::GetColor("#4477AA"), TColor::GetColor("#CCBB44"),
      TColor::GetColor("#66CCEE"), TColor::GetColor("#AA3377"),
      TColor::GetColor("#BBBBBB")};

  const int TolMutedColors[] = {
      TColor::GetColor("#88CCEE"), TColor::GetColor("#44AA99"),
      TColor::GetColor("#117733"), TColor::GetColor("#332288"),
      TColor::GetColor("#DDCC77"), TColor::GetColor("#999933"),
      TColor::GetColor("#CC6677"), TColor::GetColor("#882255"),
      TColor::GetColor("#AA4499"), TColor::GetColor("#DDDDDD")};

  const int TolLightColors[] = {
      TColor::GetColor("#BBCC33"), TColor::GetColor("#AAAA00"),
      TColor::GetColor("#77AADD"), TColor::GetColor("#EE8866"),
      TColor::GetColor("#EEDD88"), TColor::GetColor("#FFAABB"),
      TColor::GetColor("#99DDFF"), TColor::GetColor("#44BB99"),
      TColor::GetColor("#DDDDDD")};

  // From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
  const int OkabeItoColors[] = {
      TColor::GetColor("#E69F00"), TColor::GetColor("#56B4E9"),
      TColor::GetColor("#009E73"), TColor::GetColor("#F0E442"),
      TColor::GetColor("#0072B2"), TColor::GetColor("#D55E00"),
      TColor::GetColor("#CC79A7"), TColor::GetColor("#000000")};

  const int OkabeItoLightColors[] = {
      TColor::GetColor("#EBB233"), TColor::GetColor("#78C3ED"),
      TColor::GetColor("#33B18F"), TColor::GetColor("#F3E968"),
      TColor::GetColor("#338EC1"), TColor::GetColor("#DD7E33"),
      TColor::GetColor("#D694B9"), TColor::GetColor("#333333")};

  const int OkabeItoDarkColors[] = {
      TColor::GetColor("#b87f00"), TColor::GetColor("#4590ba"),
      TColor::GetColor("#007e5c"), TColor::GetColor("#c0b635"),
      TColor::GetColor("#005b8e"), TColor::GetColor("#aa4b00"),
      TColor::GetColor("#a36186"), TColor::GetColor("#333333")};

  const int  Helium_12colorScheme[] = {
      TColor::GetColor("#e25822"),  // reddishorange,
      TColor::GetColor("#DF00FF"), //'psychedelic Purple
      TColor::GetColor("#ffc922"), //'sunset yellow'
      TColor::GetColor("#87CEEB"), //'skyblue'
      TColor::GetColor("#0859C6"), //blue
      TColor::GetColor("#654522"),  // yellowishbrown,
      TColor::GetColor(43, 206, 72), //green
      TColor::GetColor("#FF0800"), // candy Apple
      TColor::GetColor("#90AD1C"),
      TColor::GetColor("#6495ED"), //'cornflowerblue'
      TColor::GetColor("#00FFFF"),//'aqua'
      TColor::GetColor("#FF0000"),//'red'
      TColor::GetColor("#32CD32"),//'limegreen'
      TColor::GetColor("#FFD700")  // Gold
      };

  const int npalettes = 13;
  static std::vector<std::vector<int> > all_palettes(npalettes);
  static bool first_call = true;

  if (first_call) {
    for (int j = 0; j < npalettes; ++j) {
      int ncolors = 0;
      const int* first_color = NULL;

      switch (j) {
        case 0:
          ncolors = 26;
          first_color = AlphabetColors;
          break;
        case 1:
          ncolors = 21;
          first_color = KellyColors;
          break;
        case 2:
          ncolors = 35;  // 36 because I commented one out
          first_color = Palette36Colors;
          break;
        case 3:
          ncolors = 30;
          first_color = GlasbeyColors;
          break;
        case 4:
          ncolors = 8;
          first_color = BrewerSet1Colors;
          break;
        case 5:
          ncolors = 8;
          first_color = BrewerDark2Colors;
          break;
        case 6:
          ncolors = 7;
          first_color = TolBrightColors;
          break;
        case 7:
          ncolors = 10;
          first_color = TolMutedColors;
          break;
        case 8:
          ncolors = 9;
          first_color = TolLightColors;
          break;
        case 9:
          ncolors = 8;
          first_color = OkabeItoColors;
          break;
        case 10:
          ncolors = 8;
          first_color = OkabeItoLightColors;
          break;
        case 11:
          ncolors = 8;
          first_color = OkabeItoDarkColors;
          break;
        case 12:
         ncolors = 14; 
         first_color = Helium_12colorScheme;
         break;
        
         default:
          assert(("palette must be 0-9", 0));
      }

      for (int i = 0; i < ncolors; ++i) {
        all_palettes[j].push_back(first_color[i]);
      }
    }
    first_call = false;
  }
  return all_palettes[palette];
}

//======================================================================
void MnvColors::AutoColorHists(TPad* pad, int palette) {
  const std::vector<int>& colors = GetColors(palette);

  std::vector<TH1*> hists = GetPadHists(pad);
  for (unsigned int i = 0; i < hists.size(); ++i) {
    hists[i]->SetLineColor(colors[i % colors.size()]);
  }
  pad->Draw();
}

//======================================================================
std::vector<TH1*> MnvColors::GetPadHists(TPad* pad) {
  std::vector<TH1*> ret;
  TIter next(pad->GetListOfPrimitives());
  TObject* obj;
  while ((obj = next())) {
    // cout << obj->GetName() << endl;
    if (obj->IsA()->InheritsFrom(TH1::Class())) {
      ret.push_back((TH1*)obj);
    }
  }
  return ret;
}

//======================================================================
const char* MnvColors::GetPaletteName(EColorPalettes e) {
  switch (e) {
    case kAlphabetPalette:
      return "Alphabet";
    case kKellyPalette:
      return "Kelly";
    case k36Palette:
      return "ThirtySix";
    case kGlasbeyPalette:
      return "Glasbey";
    case kBrewerSet1Palette:
      return "BrewerDark1";
    case kBrewerDark2Palette:
      return "BrewerDark2";
    case kTolBrightPalette:
      return "TolBright";
    case kTolMutedPalette:
      return "TolMuted";
    case kTolLightPalette:
      return "TolLight";
    case kOkabeItoPalette:
      return "OkabeIto";
    case kOkabeItoLightPalette:
      return "OkabeItoLight";
    case kOkabeItoDarkPalette:
      return "OkabeItoDark";
    case kHeliumPalette:
      return "HeliumPalette"; 
   default:
      assert(("Invalid Palette Choice.", 0));
      return "INVALID";
  }
}

#endif  // MNVCOLORS_CXX
