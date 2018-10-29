#include "Setting.h"

Setting::Setting()
{
}

//Setting::~Setting()
//{
//}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->SetTitleFont(42,"");
  h->SetTitleSize(0.04,"");

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetTH2(TH2 *h, TString name, TString xname, TString yname, double min){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->SetTitleFont(42,"");
  h->SetTitleSize(0.04,"");

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.00);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetGr(TGraphErrors *gr, TString name, TString xname, TString yname, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.00);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetGr(TGraphAsymmErrors *gr, TString name, TString xname, TString yname, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.00);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetGr(TGraph2D *gr, TString name, TString xname, TString yname, TString zname, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(1.10);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.16);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelSize(0.03);
  gr->GetYaxis()->SetLabelOffset(0.01);

  gr->GetZaxis()->SetTitle(zname);
  gr->GetZaxis()->CenterTitle();
  gr->GetZaxis()->SetTitleFont(42);
  gr->GetZaxis()->SetTitleOffset(1.00);
  gr->GetZaxis()->SetTitleSize(0.06);
  gr->GetZaxis()->SetLabelFont(42);
  gr->GetZaxis()->SetLabelSize(0.03);
  gr->GetZaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetZaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::SetTF1(TF1 *f, int LColor, int LWidth, int LStyle, int Npx){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::BinLogX(TH2 *h){
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();
   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];
   for (int i=0;i<=bins;i++){ new_bins[i] = TMath::Power(10,from + i*width); }
   axis->Set(bins, new_bins);
   delete new_bins;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::BinLogY(TH2 *h){
   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();
   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];
   for (int i=0;i<=bins;i++){ new_bins[i] = TMath::Power(10,from + i*width); }
   axis->Set(bins, new_bins);
   delete new_bins;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Setting::landaugaus(double *x, double *par) {
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 100.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val;

// MP shift correction
  mpc = par[1] - mpshift * par[0];
// Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
// Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);

     xx = xupp - (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  val = par[2] * step * sum * invsq2pi / par[3];

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Setting::Initialize(){
//  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");

  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
  gStyle->SetStatFont(42);
  gStyle->SetGridWidth(1);
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"X"); // tertiary*10000 + secondary*100 + first
  gStyle->SetNdivisions(505,"Y"); // tertiary*10000 + secondary*100 + first
  gStyle->SetNdivisions(505,"Z"); // tertiary*10000 + secondary*100 + first

// stat box
  gStyle->SetOptStat("");
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatFont(42);

// Pad
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);

// Title
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.050);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleStyle(0);

//  gStyle->SetTextFont(42);
//  gStyle->SetTextSize(0.02);
//  gStyle->SetEndErrorSize(5);
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
