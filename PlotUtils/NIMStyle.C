{
  TStyle * style = new TStyle("nim", "NIM paper plot style");

  // First, copy the stuff from the ROOT "Plain" style.
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetPadColor(0);
  style->SetCanvasColor(0);
  style->SetTitleColor(1);
  style->SetStatColor(0);

  // Set the size of the default canvas: 600x500 looks almost square.
  style->SetCanvasDefH(500);
  style->SetCanvasDefW(600);
  style->SetCanvasDefX(10);
  style->SetCanvasDefY(10);

  // Borders
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);

  // Color Scheme - Black & White
  style->SetCanvasColor(0);
  style->SetPalette(1);
  style->SetFrameBorderMode(0);
 
  // Line Widths
  style->SetFrameLineWidth(2);
  style->SetLineWidth(1);
  style->SetHistLineWidth(2);

  // Marker Styles
  style->SetMarkerStyle(20);

  // Stats
  style->SetOptStat(0000);
  style->SetOptFit(0000);
  
  // Margins
  style->SetPadTopMargin(0.15);
  style->SetPadBottomMargin(0.15);
  style->SetPadLeftMargin(0.15);
  style->SetPadRightMargin(0.15);

  // Axis label size and offset
  style->SetTitleOffset(1.5,"Y");
  style->SetTitleOffset(1.2,"X");
  style->SetTitleOffset(1.2,"Z");
  style->SetTitleSize(0.045,"XY");
  style->SetLabelSize(0.045,"XY");
  style->SetLabelFont(62);
  style->SetStripDecimals(false); // don't remove trailing zeros

  // Errors
  style->SetEndErrorSize(3);
  style->SetErrorX(0.5);

  // Legend
  style->SetLegendBorderSize(1);
  style->SetLegendFillColor(0);
  style->SetFillColor(10);
  style->SetLegendFont(62);

  // Title
  style->SetTitleBorderSize(0);
  style->SetTitleX(0.1f);
  style->SetTitleW(0.8f);
  style->SetTitleFont(62);
  style->SetTitleFontSize(0.0555555);
  style->SetOptTitle(0);

  // Finally...
  gROOT->SetStyle("nim");
  gROOT->ForceStyle();
}

