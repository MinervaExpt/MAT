// GridCanvas.cxx -
//
// Philip Rodrigues
// Thursday, February  7 2013
//

#include "GridCanvas.h"
#include "TLatex.h"
#include "TPad.h"
#include "TH1.h"
#include "TAxis.h"
#include "TList.h"
#include "TGraph.h"
#include "TStyle.h"

// from plot.h
#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "THStack.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TKey.h"

#include "TLine.h"



//#include "plot.h"

#include <iostream>

using namespace PlotUtils;

//======================================================================
// Copied from plot.h
double GridCanvas::getPadMax(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
    /*if ( obj->IsA()->InheritsFrom(THStack::Class()) ) {
      THStack* curHist=(THStack*)obj;
      const double thisMax=curHist->GetMaximum();
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }*/
  }
  return runningMax;
}

GridCanvas::GridCanvas(const char* name, int nPadsX, int nPadsY, int ww, int wh)
  : TCanvas(name, "title", ww, wh), fNPadsX(nPadsX), fNPadsY(nPadsY),
    fInterpadSpace(0.001),
    fXTitleLatex(new TLatex), fYTitleLatex(new TLatex),
    fTitleAlignment(kAlignCenter),
    fXTitle(""), fYTitle(""),
    fXTitleDrawn(false), fYTitleDrawn(false),
    fTitleFont(-1), fTitleSize(-1),
    fManualXLabels(false),
    fheadroom(1.2)
{
  fPads.resize(fNPadsX*fNPadsY);
  fPads2D.resize(fNPadsX);

  //const double padWidth=(1-fLeftMargin-fRightMargin)/fNPadsX;
  //const double padHeight=(1-fTopMargin-fBottomMargin)/fNPadsY;

  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      TPad *pad = new TPad(TString::Format("pad%d", counter), "foo",
                           0, 0, 1, 1);
      pad->SetNumber(counter+1);
      fPads[counter]=pad;
      fPads2D[i][j]=pad;

      TCanvas::cd();
      pad->Draw();
    }
  }
  ResetPads();
}

GridCanvas::~GridCanvas()
{
  delete fXTitleLatex;
  delete fYTitleLatex;
  for(unsigned int i=0; i<fPads.size(); ++i) delete fPads[i];
}

void GridCanvas::ResetPads()
{
  const double gridWidth=(1-fLeftMargin-fRightMargin);
  const double gridHeight=(1-fTopMargin-fBottomMargin);

  const double frameWidth=gridWidth/fNPadsX;
  const double frameHeight=gridHeight/fNPadsY;
  const double aspectRatio=frameWidth/frameHeight;

  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      // const double thisPadHeight = j==0 ? padHeight+fBottomMargin  : padHeight;
      // const double thisPadWidth  = i==0 ? padWidth+fLeftMargin : padWidth;

      TPad *pad = fPads[counter];

      pad->SetFillStyle(4000);

      const double left=fLeftMargin + i*frameWidth + fInterpadSpace/aspectRatio;
      const double right=fRightMargin + (fNPadsX-i-1)*frameWidth + fInterpadSpace/aspectRatio;
      const double bottom=fBottomMargin + j*frameHeight + fInterpadSpace;
      const double top=fTopMargin + (fNPadsY-j-1)*frameHeight + fInterpadSpace;

      pad->SetLeftMargin(left);
      pad->SetRightMargin(right);
      pad->SetBottomMargin(bottom);
      pad->SetTopMargin(top);

      //printf("(%d, %d) c=%d l=%.2f r=%.2f b=%.2f t=%.2f\n", i, j, counter, left, right, bottom, top);
      // pad->SetBottomMargin(j==0 ? (thisPadHeight-padHeight)/thisPadHeight : fInterpadSpace);
      // pad->SetTopMargin(fInterpadSpace);
      // pad->SetLeftMargin(i==0 ? (thisPadWidth-padWidth)/thisPadWidth : fInterpadSpace);
      // pad->SetRightMargin(fInterpadSpace);
    }
  }
}

TH1* GridCanvas::GetPadHist(TPad* pad)
{
  TList* prims=pad->GetListOfPrimitives();
  for(int i=0; i<prims->GetEntries(); ++i){
    TObject* obj=prims->At(i);
    // std::cout << "Object name=" << obj->GetName() << " classname=" << obj->ClassName() << std::endl;
    TH1* h=dynamic_cast<TH1*>(obj);
    if(h) return h;
    TGraph* gr=dynamic_cast<TGraph*>(obj);
    if(gr) return gr->GetHistogram();
  }
  return 0;
}

/*THStack* GridCanvas::GetPadStack(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(THStack::Class()) ) {
      THStack* curHist=(THStack*)obj;
      return curHist;
    }
  }
  return NULL;
}*/

void GridCanvas::Paint(Option_t* option)
{
  bool anyPadModified=false;
  for(unsigned int i=0; i<fPads.size(); ++i) anyPadModified=anyPadModified || fPads[i]->IsModified();
  if(anyPadModified || IsModified()){
    SetHistTexts();
  }
  TCanvas::Paint(option);
}

void GridCanvas::SetHistTexts()
{
  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      // const double thisPadHeight = j==0 ? padHeight+fBottomMargin  : padHeight;
      // const double thisPadWidth  = i==0 ? padWidth+fLeftMargin : padWidth;

      TPad *pad = fPads[counter];

      TH1* hist=GetPadHist(pad);
      if(!hist) continue;

      hist->GetXaxis()->SetTitleSize(0);
      hist->GetYaxis()->SetTitleSize(0);

      if(i!=0) hist->GetYaxis()->SetLabelSize(0);
      if(j!=0 || fManualXLabels) hist->GetXaxis()->SetLabelSize(0);
    }
  }
  if(fXTitle=="" && GetPadHist(fPads[0])) SetXTitle(GetPadHist(fPads[0])->GetXaxis()->GetTitle());
  if(fYTitle=="" && GetPadHist(fPads[0])) SetYTitle(GetPadHist(fPads[0])->GetYaxis()->GetTitle());

  DrawTitles();
}

void GridCanvas::SetManualXLabels(int nLabels, const double* positions, const char** valueStrings,
                                  double yoffset)
{
  fManualXLabels=true;
  for(int i=0; i<fNPadsX; ++i){
    std::cout << "PAD " << i << std::endl;
    fPads2D[i].resize(fNPadsY);
    int j=0;
    int counter=fNPadsX*(fNPadsY-1-j)+i;
    
    TPad *pad = fPads[counter];
    pad->cd();
    pad->Update();

    TH1* hist=GetPadHist(pad);
    if(hist==NULL) continue;
    double x1=pad->GetUxmin();
    double x2=pad->GetUxmax();
    double y1=pad->GetUymin();
    double y2=pad->GetUymax();

    double lmarg=pad->GetLeftMargin();
    double rmarg=pad->GetRightMargin();
    double tmarg=pad->GetTopMargin();
    double bmarg=pad->GetBottomMargin();

    // std::cout << "y2=" << y2 << " y1=" << y1 << " ypos=" << (y1-yoffset*(y2-y1)) << std::endl;

    for(int i=0; i<nLabels; ++i){
      // We have to place the TLatex in NDC so it doesn't move when
      // the user changes the y axis limits. This'll still break if
      // the user changes the x limits, but you can't have
      // everything...
      // 
      // Presumably the real fix is to paint() the latex in
      // GridCanvas::Paint() instead of doing this, but laziness

      double y=y1-yoffset*(y2-y1);

      double xndc= lmarg + (positions[i]-x1)*(1-rmarg-lmarg)/(x2-x1);
      double yndc= bmarg + (y-y1)*(tmarg-bmarg)/(y2-y1);

      TLatex* la=new TLatex(xndc, yndc, /* positions[i], y1-yoffset*(y2-y1), */
                            valueStrings[i]);
      la->SetNDC();
      la->SetTextAlign(23);
      
      la->SetTextFont(hist->GetXaxis()->GetLabelFont());
      // Can't do this because we set the label size to zero so ROOT's
      // labels aren't shown
      // la->SetTextSize(axis->GetLabelSize());
      la->SetTextSize(gStyle->GetLabelSize());
      la->Draw();
    }
    std::cout << "DONE" << std::endl;
  }
}

void GridCanvas::SetXTitle(const char* xtitle)
{
  fXTitle=xtitle;
}

void GridCanvas::SetYTitle(const char* ytitle)
{
  fYTitle=ytitle;
}

TLatex* GridCanvas::GetXTitle()
{
  return fXTitleLatex;
}

TLatex* GridCanvas::GetYTitle()
{
  return fYTitleLatex;
}

void GridCanvas::SetXLimits(double xmin, double xmax)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    if(!h) continue;
    h->GetXaxis()->SetRangeUser(xmin, xmax);
  }
}

void GridCanvas::SetYLimits(double ymin, double ymax)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    if(!h) 
    {
      /*THStack* hs=GetPadStack(fPads[i]);
      std::cout<<i<<", "<<hs<<std::endl;
      if (hs)
      {
        hs->SetMinimum(ymin);
        hs->SetMaximum(ymax);
      }*/
      // This is fine except for the case where there's a THStack
      continue;
    }
    h->GetYaxis()->SetRangeUser(ymin, ymax);
  }
}

void GridCanvas::DrawTitles()
{
  fXTitleLatex->SetTitle(fXTitle);
  fYTitleLatex->SetTitle(fYTitle);

  if(GetPadHist(fPads[0])){
    // TH1* h=GetPadHist(fPads[0]);
    // fXTitleLatex->SetTextFont(h->GetXaxis()->GetTitleFont());
    // fYTitleLatex->SetTextFont(h->GetYaxis()->GetTitleFont());

    // fXTitleLatex->SetTextSize(h->GetXaxis()->GetTitleSize());
    // fYTitleLatex->SetTextSize(h->GetYaxis()->GetTitleSize());

    fXTitleLatex->SetTextFont(fTitleFont==-1 ? 43 : fTitleFont );
    fYTitleLatex->SetTextFont(fTitleFont==-1 ? 43 : fTitleFont);

    fXTitleLatex->SetTextSize(fTitleSize==-1 ? 30 : fTitleSize);
    fYTitleLatex->SetTextSize(fTitleSize==-1 ? 30 : fTitleSize);
  }

  fYTitleLatex->SetTextAngle(90);

  double xposx, xposy, yposx, yposy;
  switch(fTitleAlignment){
  case kAlignRight:
    xposx=1-fRightMargin;
    yposx=0.02;

    xposy=0.02;
    yposy=1-fTopMargin;

    fXTitleLatex->SetTextAlign(31);
    fYTitleLatex->SetTextAlign(33);
    break;
  case kAlignCenter:
    xposx=fLeftMargin+0.5*(1-fLeftMargin-fRightMargin);
    yposx=0.02;

    xposy=0.01;
    yposy=fBottomMargin+0.5*(1-fLeftMargin-fRightMargin);

    fXTitleLatex->SetTextAlign(21);
    fYTitleLatex->SetTextAlign(23);
    break;
  default:
    std::cerr << "Unknown alignment type in GridCanvas::DrawTitles()" << std::endl;
    return;
  }
  fXTitleLatex->SetX(xposx);   fXTitleLatex->SetY(yposx);
  fYTitleLatex->SetX(xposy);   fYTitleLatex->SetY(yposy);

  fXTitleLatex->SetNDC();
  fYTitleLatex->SetNDC();

  if(!fXTitleDrawn){
    std::cout << "Drawing x title" << std::endl;
    TCanvas::cd();
    fXTitleLatex->Draw();
    fXTitleDrawn=true;
  }

  if(!fYTitleDrawn){
    std::cout << "Drawing y title" << std::endl;
    TCanvas::cd();
    fYTitleLatex->Draw();
    fYTitleDrawn=true;
  }
}

void GridCanvas::SetLogx(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogx(value);
}

void GridCanvas::SetLogy(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogy(value);
}

void GridCanvas::SetLogz(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogz(value);
}

void GridCanvas::SetGridx(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetGridx(value);
}

void GridCanvas::SetGridy(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetGridy(value);
}

void GridCanvas::SetTicksy(const char* option)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    h->GetYaxis()->SetTicks(option);
    if(strcmp(option, "+")==0){
      // TODO: Work out what the hell this number means, why it's so
      // stupid, and how to make it leave the numbers unchanged
      h->GetYaxis()->SetLabelOffset(-0.015);
    }
  }
}

void GridCanvas::SetTicksx(const char* option)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    h->GetXaxis()->SetTicks(option);
    if(strcmp(option, "+")==0){
      // TODO: Work out what the hell this number means, why it's so
      // stupid, and how to make it leave the numbers unchanged
      h->GetXaxis()->SetLabelOffset(-0.015);
    }
  }
}

void GridCanvas::Remax(double ymin)
{
  double allMax=-9e99;
  for(unsigned int i=0; i<fPads.size(); ++i) allMax=std::max(allMax, getPadMax(fPads[i]));
  SetYLimits(ymin, fheadroom*allMax);
}

void GridCanvas::SetLeftMargin(Float_t margin)
{
  TCanvas::SetLeftMargin(margin);
  ResetPads();
}

void GridCanvas::SetRightMargin(Float_t margin)
{
  TCanvas::SetRightMargin(margin);
  ResetPads();
}

void GridCanvas::SetHeadroom(double headroom)
{
  fheadroom = headroom;
  ResetPads();
}

void GridCanvas::SetTopMargin(Float_t margin)
{
  TCanvas::SetTopMargin(margin);
  ResetPads();
}

void GridCanvas::SetBottomMargin(Float_t margin)
{
  TCanvas::SetBottomMargin(margin);
  ResetPads();
}

// Return a unique string on every call, so we can name root hists
// without it clobbering them
TString GridCanvas::uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}


TH1* GridCanvas::GetProjection(TH2D* h2d, std::vector<int> bins, int index, bool ydirection, bool widthscale, double* multipliers)
{
    TH1* proj;
    int bin_max = bins[index+1] - 1;
    int bin_min = bins[index];
    if (ydirection) proj = h2d->ProjectionY(uniq(), bin_min, bin_max);
    else proj = h2d->ProjectionX(uniq(), bin_min, bin_max);
    proj->SetLineStyle(h2d->GetLineStyle());
    proj->SetLineWidth(h2d->GetLineWidth());
    proj->SetMarkerStyle(h2d->GetMarkerStyle());
    proj->SetMarkerSize(h2d->GetMarkerSize());
    proj->GetXaxis()->SetNdivisions(4);
    proj->GetYaxis()->SetNdivisions(4);
        
    if (widthscale) proj->Scale(1, "width");
    if (multipliers) proj->Scale(multipliers[index]);
    
    return proj;
}

void GridCanvas::DrawOneHist(TH2D* h2d, const char* opts, bool widthscale, std::vector<int> bins, bool ydirection, double* multipliers)
{
  for (unsigned int i = 0; i < bins.size() - 1; i++)
  {
    cd(i+1);
    auto proj = GetProjection(h2d, bins, i, ydirection, widthscale, multipliers);
    proj->Draw(opts);
  }
  Remax();
}

void GridCanvas::DrawRatio(TH2D* numerator, TH2D* denominator, const char* opts, std::vector<int> bins, bool ydirection, bool include_denom_error)
{
  // IMPORTANT: If drawing data/MC for a xsec, set include_denom_error to false;
  for (unsigned int i = 0; i < bins.size() - 1; i++)
  {
    cd(i+1);
    auto proj1 = GetProjection(numerator, bins, i, ydirection, false, NULL);
    auto proj2 = GetProjection(denominator, bins, i, ydirection, false, NULL);
    if (!include_denom_error)
    {
      for (int bin = 0; bin < proj2->GetNbinsX(); bin++)
        proj2->SetBinError(bin, 0);
    }
    proj1->Divide(proj2); 
    TLine *line = new TLine(proj1->GetBinLowEdge(1),1,proj1->GetBinLowEdge(proj1->GetNbinsX() + 1),1);
    line->SetLineStyle(2);
    proj1->Draw(opts); 
    line->Draw();
    proj1->Draw(opts); 
  }
  Remax();
}

void GridCanvas::DrawStack(std::vector<TH2D*> hists, const char* opts, bool widthscale, std::vector<int> bins, bool ydirection, double* multipliers)
{
  for (unsigned int i = 0; i < bins.size() - 1; i++)
  {
    cd(i+1);
    //THStack* stack = new THStack();
    std::vector<TH1*> projections;
    for (int j = hists.size() - 1; j >= 0; j--)
    {
      auto proj = GetProjection(hists[j], bins, i, ydirection, widthscale, multipliers);
      for (unsigned int I = 0; I < projections.size(); I++) projections[I]->Add(proj);
      projections.push_back(proj);
    }
    
    for (unsigned int j = 0; j < hists.size(); j++)
    {
      TString newopts(opts);
      if (j != 0) newopts += " same";
      projections[j]->Draw(newopts);
    }
  }
  Remax();
}

void GridCanvas::DrawMultipliers(int n_multipliers, double* multipliers)
{
    for (int i = 0; i < n_multipliers; i++)
    {
        if(multipliers && multipliers[i]!=1)
        {
            auto pad = cd(i+1);
            TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                                 1-pad->GetTopMargin()-0.08,
                                 TString::Format("#times %.1f", multipliers[i]));
            la2->SetTextAlign(33); // top right
            la2->SetNDC();
            la2->SetTextFont(42);
            la2->SetTextSize(0.035);
            la2->Draw();
        }
    }
}

void GridCanvas::DrawBinRanges(TH2* h, int axis, std::vector<int> bins, const char* varName, double textsize, const char* numFormatStr, int position)
{
    for (unsigned int i = 0; i < bins.size() - 1; i++)
    {
        cd(i+1);
        int bin_max = bins[i+1] - 1;
        int bin_min = bins[i];
        DrawBinRange(h, axis, bin_min, bin_max, varName, textsize, numFormatStr, position);
    }
}

void GridCanvas::DrawBinRange(TH2* h, int axis, int bin_min, int bin_max, const char* varName, double textsize, const char* numFormatStr, int position)
{

  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin_min) : h->GetYaxis()->GetBinLowEdge(bin_min);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin_max) :  h->GetYaxis()->GetBinUpEdge(bin_max);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));

  if(position == 3){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else if (position == 2)
  {
    la=new TLatex((1-gPad->GetRightMargin()-0.02 + gPad->GetLeftMargin()+0.02) * 0.5,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(23); // top center
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  //Okay math time. The grid will be 150 * ncolumns wide and 150*nrows tall
  //SetTextSize is based off the aspect ratio and then the character height is defined as (size)*pad_width (tall)or (size)*pad_height (wide)
  // This is a problem because the pads for all "frames" are the full canvas with modified margins. That means the size of the text stays constant for any number of "frames". We want smaller text for more frames.
  // A good size for wide plots has been 24. 800 pixel wide * 0.03
  // So, 150*ncolumn = width we want 24 -> 24/(150*ncolumn) or 24/(150*nrow)
  /*if(gridx==0||gridy==1) la->SetTextSize(0.03);
  else if(gridx<gridy) la->SetTextSize(28/(150.0*gridy));//Tall
  else if(gridx==gridy) la->SetTextSize(10/(150.0*gridx));//Wide
  else la->SetTextSize(18/(150.0*gridx));//Wide*/
  la->SetTextSize(textsize);
  la->Draw();
}

// https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html
//ClassImp(GridCanvas)
