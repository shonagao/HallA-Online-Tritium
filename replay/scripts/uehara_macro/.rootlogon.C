{
  gROOT->SetStyle("Plain");
  gROOT->LoadMacro("/adaqfs/home/a-onl/tritium_work/uehara/rootmacro/online_check/.rootmacro.cc");
  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
  gStyle->SetTextFont(42);
  gStyle->SetGridWidth(0);
  gStyle->SetNdivisions(505); // tertiary*10000 + secondary*100 + first

// Stat box
  gStyle->SetOptStat("ei");
  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFont(42);
// Pad
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);
// Title
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleFont(42,"");
// Label
  gStyle->SetStripDecimals(kFALSE);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.012,"X");
  gStyle->SetLegendFont(42);

  const Int_t NRGBs = 5;
  const Int_t NCont = 99;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}
