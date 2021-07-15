#ifndef GRIDCANVAS_H
#define GRIDCANVAS_H

#include "TCanvas.h"
#include "TString.h"

#include <vector>

class TLatex;
class TPad;
class TH1;
class TH2;
class TH2D;
class THStack;

namespace MAT{
class GridCanvas : public TCanvas
{
public:
  GridCanvas() {}; // To shut ROOT up
  GridCanvas(const char* name, int nPadsX, int nPadsY, int ww=700, int wh=500);
  
  
  virtual ~GridCanvas();

  void SetXTitle(const char* xtitle);
  void SetYTitle(const char* ytitle);

  TLatex* GetXTitle();
  TLatex* GetYTitle();

  std::vector<TPad*> GetPads() const { return fPads; }


  void SetHistTexts();

  void SetXLimits(double xmin, double xmax);
  void SetYLimits(double ymin, double ymax);

  enum ETitleAlignment { kAlignRight, kAlignCenter };

  void SetTitleAlignment(ETitleAlignment alignment) { fTitleAlignment=alignment; }

  void SetInterpadSpace(double space) { fInterpadSpace=space; ResetPads(); }

  void Paint(Option_t*);

  virtual void SetLogx(Int_t value=1);
  virtual void SetLogy(Int_t value=1);
  virtual void SetLogz(Int_t value=1);

  virtual void SetGridx(Int_t value=1);
  virtual void SetGridy(Int_t value=1);

  virtual void SetLeftMargin(Float_t margin);
  virtual void SetRightMargin(Float_t margin);
  virtual void SetTopMargin(Float_t margin);
  virtual void SetBottomMargin(Float_t margin);

  void SetTicksx(const char* option);
  void SetTicksy(const char* option);

  void ResetPads();

  void SetTitleFont(int font) { fTitleFont=font; }
  int  GetTitleFont() const { return fTitleFont; }

  void   SetTitleSize(double size) { fTitleSize=size; }
  double GetTitleSize() const { return fTitleSize; }

  void Remax(double ymin=0);

  void SetManualXLabels(int nLabels, const double* positions, const char** valueStrings, double yoffset=0.1);
  
  void DrawOneHist(TH2D* hist, const char* opts, bool widthscale, std::vector<int> bins, bool direction, double* multipliers=NULL);
  void DrawStack(std::vector<TH2D*> hists, const char* opts, bool widthscale, std::vector<int> bins, bool direction, double* multipliers=NULL);
  void DrawRatio(TH2D* numerator, TH2D* denominator, const char* opts, std::vector<int> bins, bool ydirection, bool include_denom_error);
  
  void DrawMultipliers(int n_multipliers, double* multipliers);
  
  void DrawBinRanges(TH2* h, int axis, std::vector<int> bins, const char* varName, double textsize = 0.1, const char* numFormatStr=".2f", int position = 2);
  
  void SetHeadroom(double headroom);
  

private:

  TH1* GetProjection(TH2D* hist, std::vector<int> bins, int index, bool direction, bool widthscale, double* multipliers=NULL);
  
  void DrawBinRange(TH2* h, int axis, int bin_min, int bin_max, const char* varName, double textsize = 0.1, const char* numFormatStr=".2f", int position = 2);

  TString uniq();

  double getPadMax(TPad* pad);

  TH1* GetPadHist(TPad* pad);
  
  THStack* GetPadStack(TPad* pad);

  void DrawTitles();

  int fNPadsX, fNPadsY;
  std::vector<TPad*> fPads;
  std::vector< std::vector<TPad*> > fPads2D;
  double fInterpadSpace;
  TLatex* fXTitleLatex;
  TLatex* fYTitleLatex;

  ETitleAlignment fTitleAlignment;
  TString fXTitle, fYTitle;
  bool fXTitleDrawn, fYTitleDrawn;
  
  int fTitleFont;
  double fTitleSize;

  bool fManualXLabels;
  
  double fheadroom;

  // https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html
  ClassDef(GridCanvas, 0);
};
}

#endif
