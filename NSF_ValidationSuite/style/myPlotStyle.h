#ifndef MYPLOTSTYLE_H
#define MYPLOTSTYLE_H

void myPlotStyle()
{
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);

  // use large fonts

  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);

  gStyle->SetTextSize(0.08);

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetTitleFontSize(0.08);
  //gStyle->SetTitleOffset(1.0, "Y");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.1, "Z");
  //gStyle->SetLabelOffset(0.02, "X");
  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.02, "Y");
  gStyle->SetTitleAlign(23);

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTickLength(0.01, "Y");
  gStyle->SetTickLength(0.02, "X");

  gStyle->SetNdivisions(505, "XYZ");
  gStyle->SetStripDecimals(false);

  gStyle->SetPalette(54); // kBlueYellow

  TGaxis::SetMaxDigits(2);
}

#endif
