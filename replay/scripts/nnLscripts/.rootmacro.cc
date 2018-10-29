const double PI = 4.*atan(1.);
const double deg_to_rad = PI / 180.;
const double rad_to_deg = 180. / PI;
const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));
const double fwhm_to_sigma = 1./sigma_to_fwhm;
const double cm_to_barn = 1e+24;
const double alpha = 1./137.035999074; // fine structure constant
const double hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
const double hbarc = 197.3269718;      // conversion constant (MeV x fm)
const double kb = 8.6173324*1e-5;      // Boltzmann constant
const double e = 1.602176565*1e-19;    // electron charge magnitude (C)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const double re = 2.817e-13;           // classical electron radius (cm)
const double Na = 6.02214129*1e+23;    // Avogadro constant
const double Me = 0.510998928;         // electron     mass (MeV/c2)
const double Mmu = 105.6583715;        // muon         mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double Mpi0 = 134.9766;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mn = 939.565379;          // proton       mass (MeV/c2)
const double Mu = 931.494061;          // proton       mass (MeV/c2)
const double ML = 1115.683;            // Lambda       mass (MeV/c2)
const double MS0 = 1192.642;           // Sigma Zero   mass (MeV/c2)
const double MSm = 1197.449;           // Sigma Minus  mass (MeV/c2)
const double MSp = 1189.37;            // Sigma Plus   mass (MeV/c2)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor=1, int FStyle=0, int FColor=0){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetTitleSize(0.04,"");
  h->SetTitleFont(42,"");
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetLabelSize(0.045);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(1);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetLabelSize(0.045);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);
  gr->GetXaxis()->SetLabelSize(0.045);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.20);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraphErrors *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);
  gr->GetXaxis()->SetLabelSize(0.045);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.20);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraphAsymmErrors *gr, TString name, TString xname, TString yname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);
  gr->GetXaxis()->SetLabelSize(0.045);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.20);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraph2D *gr, TString name, TString xname, TString yname, TString zname, int LColor=1, int LWidth=1, int LStyle=1, int MColor=1, int MStyle=20, double MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(1.10);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);
  gr->GetXaxis()->SetLabelSize(0.045);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.16);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetLabelSize(0.045);

  gr->GetZaxis()->SetTitle(zname);
  gr->GetZaxis()->CenterTitle();
  gr->GetZaxis()->SetTitleFont(42);
  gr->GetZaxis()->SetTitleOffset(1.00);
  gr->GetZaxis()->SetTitleSize(0.06);
  gr->GetZaxis()->SetLabelFont(42);
  gr->GetZaxis()->SetLabelOffset(0.01);
  gr->GetZaxis()->SetLabelSize(0.045);
  ((TGaxis*)gr->GetZaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetLegend(TLegend *l, int FStyle=1001, int FColor=0, double TSize=0.050, int NColumn=1){
  l->SetFillStyle(FStyle);
  l->SetFillColor(FColor);
  l->SetTextSize(TSize);
  l->SetNColumns(NColumn);
}


