//enum EColor { kWhite =0,   kBlack =1,   kGray=920,
//              kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
//              kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

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
  h->SetLineWidth(0);
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

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.00);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
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

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
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

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.16);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);

  gr->GetZaxis()->SetTitle(zname);
  gr->GetZaxis()->CenterTitle();
  gr->GetZaxis()->SetTitleFont(42);
  gr->GetZaxis()->SetTitleOffset(1.00);
  gr->GetZaxis()->SetTitleSize(0.06);
  gr->GetZaxis()->SetLabelFont(42);
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
void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// void BinLogX(TH2 *h){
//    TAxis *axis = h->GetXaxis(); 
//    int bins = axis->GetNbins();
//    Axis_t from = axis->GetXmin();
//    Axis_t to = axis->GetXmax();
//    Axis_t width = (to - from) / bins;
//    Axis_t *new_bins = new Axis_t[bins + 1];
//    for (int i=0;i<=bins;i++){ new_bins[i] = TMath::Power(10,from + i*width); }
//    axis->Set(bins, new_bins);
//    delete new_bins; 
// }  

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// void BinLogY(TH2 *h){
//    TAxis *axis = h->GetYaxis(); 
//    int bins = axis->GetNbins();
//    Axis_t from = axis->GetXmin();
//    Axis_t to = axis->GetXmax();
//    Axis_t width = (to - from) / bins;
//    Axis_t *new_bins = new Axis_t[bins + 1];
//    for (int i=0;i<=bins;i++){ new_bins[i] = TMath::Power(10,from + i*width); }
//    axis->Set(bins, new_bins);
//    delete new_bins;
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double landaugaus(double *x, double *par) {
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double landaugaus_pol0(double *x, double *par) {
  return landaugaus(x,par) + par[4];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double landaugaus_pol1(double *x, double *par) {
  return landaugaus(x,par) + par[4] + par[5]*x[0];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double landaugaus_pol2(double *x, double *par) {
  return landaugaus(x,par) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
//  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + par[4]*TMath::Gaus(x[0],par[3],par[2],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expexpgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + expgaus(x,&par[4]) + par[8]*TMath::Gaus(x[0],par[3],par[2],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//struct Vavilov_Func {
//   Vavilov_Func() {}
//
//   double operator() (const double *x, const double *p) {
//      double kappa = p[0];
//      double beta2 = p[1];
//      return p[4]*( pdf.Pdf( (x[0]-p[2])/p[3], kappa,beta2) );
//   }
//
//   ROOT::Math::VavilovAccurate pdf;
//};
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//double Vavilov(double *x, double *par){
//  ROOT::Math::VavilovAccurate pdf;
////  double val = par[2]*TMath::Vavilov((x[0]-par[3]),par[0],par[1]);
//  double val = par[4]*( pdf.Pdf((x[0]-par[2])/par[3],par[0],par[1]));
//  return val;
//}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Fermi(double *x, double *par){
  // par[0]=Cliff Pos(x[0]/par[0]==1)) (should be > 0)
  // par[1]=Shape parameter
  // par[2]=Scale factor (should be > 0)
  double val = par[2] * ( 1. / ( exp((x[0]-par[0])/par[1]) + 1. ) );
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Pois( double *x, double *par )
{
//  par[0] : lambda, average of #photon
//  par[1] : energy resolution factor
//  par[2] : menseki
  double val = 0;
  for( int np=1; np<500; np++ ){
    double pois;
    double sigma = par[1]*sqrt(np);
    if(np<50){
      pois = par[2]*pow( par[0], np )*exp(-par[0])/TMath::Gamma(np+1);
    }
    else{ // stirling's approximation
      pois = par[2]*pow( par[0]/np, np )*exp(-par[0]+np)/(sqrt(2.*PI)*pow(np,0.5));
    }
    val += pois/(sqrt(2.*PI)*sigma)*exp( -pow(x[0]-np,2)/2./sigma/sigma ); // adding gaussian distribution
  }
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Voigt( double *x, double *par )
{
  // par[0] : area
  // par[1] : location
  // par[2] : gaussian sigma
  // par[3] : lorentz fwhm
  double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetDivError( double x1, double ex1, double x2, double ex2 ) // x1/x2
{
  return sqrt( (1./x2*ex1)*(1./x2*ex1) + (x1/x2/x2*ex2)*(x1/x2/x2*ex2) );
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetMulError( double x1, double ex1, double x2, double ex2 ) // x1*x2
{
  return sqrt( (x2*ex1)*(x2*ex1) + (x1*ex2)*(x1*ex2) );
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetRMSError( double x1, double ex1, double x2, double ex2 ) // sqrt(x1*x1+x2*x2)
{
  return sqrt( ((x1*ex1)*(x1*ex1) + (x2*ex2)*(x2*ex2)) / (x1*x1 + x2*x2) );
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetBinomialError( int n, double p ) // statistics n, probability p
{
  return sqrt( 1./n * p * (1.-p) );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetdEdx_Scinti( double p, double M )
{
  double cons = 2.*PI*Na*re*re*Me;
  double ro = 1.032;              // density (g/cm3)
  double thick = 1.0;  // thickness (cm)
  double Z = 6.5;                 // atomic number of absorbing material
  double A = 13;               // atomic weight of absorbing material
  double z = 1.;                // charge of incident particle
  double s = Me/M;
  double C0, a, m, X0, X1;               // parameters of delta
  double Ne = Na * ro * Z / A;           // density of electrons
  double mup = sqrt(80.617 * 1e+6 * Ne); // plasma frequency (cm^3/2 x Hz)
  double delta, W, C, I, nonrela, rela;
  double dedx_wo, dedx, Delta_b, kappa;

  I = 64.7 * 1e-6;  C0 = -3.20;
  a = 0.1610;  m = 3.24;
  X1 = 2.49;  X0 = 0.1464;

  double beta = p / sqrt(p*p + M*M);
  double gamma = sqrt(p*p + M*M) / M;
  double eta = p / M;
  double X = log10(beta*gamma);

  W = 2.*Me*eta*eta / (1.+2.*s*sqrt(1+eta*eta)+s*s);
  nonrela = cons * ro * Z/A * z*z/beta/beta;
  rela = log(2.*Me*gamma*gamma*beta*beta*W/I/I)-2.*beta*beta;
  if(X<X0) delta = 0;
  else if(X>=X0 && X<X1) delta = 4.6052*X + C0 + a*pow(X1-X,m);
  else delta = 4.6052*X + C0;
  if(eta<0.1) C = 0;
  else C = (0.422377*pow(eta,-2) + 0.0304043*pow(eta,-4) - 0.00038106*pow(eta,-6)) * 1e-6 * pow(I*1e+6,2)
          +(3.850190*pow(eta,-2) - 0.1667989*pow(eta,-4) + 0.00159755*pow(eta,-6)) * 1e-9 * pow(I*1e+6,3);
  dedx_wo = nonrela * rela;
  dedx    = nonrela * (rela - delta - 2.*C/Z);
  Delta_b = nonrela * thick;
  kappa   = Delta_b / W;
// thin target, ignored logalism term
//  double  Xi = Delta_b * (log(Delta_b) - LnEps + 0.198 - delta);

  return dedx;

}


