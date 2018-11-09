/////////////////////////////
//cointime analysis        //
// by Y. Toyama Nov. 2018  //
/////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom.h"

#include "Tree.h"
#include "Setting.h"
#include "ParamMan.h"
#include "define.h"

#define Calibration
#define ALL_TRACK
#undef ALL_TRACK

const int NCanvas = 6;//num of canvas
double coffset  = 150.41;//Fbus
double ccoffset = 146.71;//Fbus
double f1offset  =  6.51057;//F1TDC
double f1coffset =  3.37;//F1TDC

//ac1 study (pion selection by cointime cut)
const int nhist_a1_asump_pi = 40;// num of a1_asump_pi histgrams projected from x vs a1_asum_pi
double nbin_a1_asump_pi = 400 / nhist_a1_asump_pi;//num of bins projected together to make a1_asum_pi histgrams

class ana_cointime : public Tree
{
 public:
         ana_cointime();
        ~ana_cointime();
  void makehist();
  void loop();
  void fit();
  void projection();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void SetOutputRoot(string ofname){output_root = ofname;}
  void SetRunList(string ifname);
  void SetInputParam(string ifname);
  Setting *set;
  ParamMan *param;

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;

    //General
    //F1
    TH1F *h_s2coin_f1,  *h_s0coin_f1; 
    TH1F *h_s2ccoin_f1, *h_s0ccoin_f1;
    TH1F *h_s2ccoin_f1ac;
    TH1F *h_s2ccoin_f1acz;
    TH1F *h_s2ccoin_f1z;
    TH1F *h_s2ccoin_f1ac_rebin;
    TH2F *h2_s2coin_f1_lpath, *h2_s0coin_f1_lpath; 
    TH2F *h2_s2coin_f1_rpath, *h2_s0coin_f1_rpath; 
    TH2F *h2_s2ccoin_f1_a1, *h2_s2ccoin_f1_a2;
    TH2F *h2_s2ccoin_f1_a1_wa2cut, *h2_s2ccoin_f1_a2_wa1cut;
    TH2F *h2_s2ccoin_f1_lpath, *h2_s0ccoin_f1_lpath; 
    TH2F *h2_s2ccoin_f1_rpath, *h2_s0ccoin_f1_rpath; 
    //Fbus
    TH1F *h_s2coin_fb, *h_s0coin_fb;
    TH1F *h_s2ccoin_fb, *h_s0ccoin_fb; 
    TH1F *h_s2ccoin_fbac;
    TH1F *h_s2ccoin_fbacz;
    TH1F *h_s2ccoin_fbac_rebin;
    TH2F *h2_s2coin_fb_lpath, *h2_s0coin_fb_lpath;
    TH2F *h2_s2coin_fb_rpath, *h2_s0coin_fb_rpath;
    TH2F *h2_s2ccoin_fb_a1, *h2_s2ccoin_fb_a2;
    TH2F *h2_s2ccoin_fb_lpath, *h2_s0ccoin_fb_lpath;
    TH2F *h2_s2ccoin_fb_rpath, *h2_s0ccoin_fb_rpath;

    TH1F *h_tof, *h_beta, *h_msq, *h_mom, *h_path, *h_t_at_targ;
    TH2F *h2_tof_beta, *h2_beta, *h2_msq_beta, *h2_mom_beta, *h2_path_beta, *h2_rf_path, *h2_rf_s2p;

    TH1F *h_Lz, *h_Rz;
    TH2F *h2_LRz;

    //Each S2 paddle
    TH1F *h_s2r_tof[RS2], *h_s2r_beta[RS2], *h_s2r_msq[RS2];
    TH1F *h_s2segcoin_f1l[16], *h_s2segcoin_fbl[16],  *h_s2segcoin_f1r[16], *h_s2segcoin_fbr[16];
    TH2F *h2_s2r_tof_beta[RS2], *h2_s2r_msq_beta[RS2], *h2_s2r_mom_beta[RS2], *h2_s2r_path_beta[RS2];
    TH2F *h2_s2segcoin_path_f1l[16], *h2_s2segcoin_path_fbl[16],  *h2_s2segcoin_path_f1r[16], *h2_s2segcoin_path_fbr[16];

    int run_num;
    TCanvas *c[NCanvas];

    TFile *ofp;
    string output_root;
    TF1 *ga_beta[16];
    TF1 *ga_coin;
    TGraphErrors *tg_beta_pos, *tg_beta_wid;
    TH2F *h_frame[2];

    double beta_pos[16],beta_wid[16];
    double ebeta_pos[16],ebeta_wid[16];

    //ac1 study (pion selection by cointime cut)
    TH2F *h2_a1_asump_x_pi, *h2_a1_asump_y_pi;
    TH1D *h_a1_asump_pi[nhist_a1_asump_pi];
    TGraph *g_a1_eff_pi;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ana_cointime::ana_cointime()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMenr");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,800 );
  }

  set = new Setting();
}
////////////////////////////////////////////////////////////////////////////
ana_cointime::~ana_cointime(){
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}
/////////////////////////////
void ana_cointime::SetRunList(string ifname){
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
    cout<<buf<<endl;
  }

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::SetInputParam(string ifname){
  param = new ParamMan(ifname.c_str());
  if(param -> SetVal())cout<<"F1TDC parameter setted"<<endl; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ana_cointime::makehist(){
  ofp = new TFile(output_root.c_str(),"recreate");
  h_s2coin_f1    = new TH1F("h_s2coin_f1"    ,"h_s2coin_f1"    ,1000,  -15, 15);
  h_s0coin_f1    = new TH1F("h_s0coin_f1"    ,"h_s0coin_f1"    ,1000,-1000,1000);
  h_s2ccoin_f1   = new TH1F("h_s2ccoin_f1"   ,"h_s2ccoin_f1"   ,1000,  -15,  15);
  h_s2ccoin_f1ac = new TH1F("h_s2ccoin_f1ac" ,"h_s2ccoin_f1ac" ,1000,  -15,  15);
  h_s2ccoin_f1acz= new TH1F("h_s2ccoin_f1acz","h_s2ccoin_f1acz",1000,  -15,  15);
  h_s2ccoin_f1z  = new TH1F("h_s2ccoin_f1z"  ,"h_s2ccoin_f1z"  ,1000,  -15,  15);
  h_s0ccoin_f1   = new TH1F("h_s0ccoin_f1"   ,"h_s0ccoin_f1"   ,1000,-10,1000);
  h_s2coin_fb    = new TH1F("h_s2coin_fb"    ,"h_s2coin_fb"    , 600,  -15,  15);
  h_s0coin_fb    = new TH1F("h_s0coin_fb"    ,"h_s0coin_fb"    ,1000,-1000,1000);
  h_s2ccoin_fb   = new TH1F("h_s2ccoin_fb"   ,"h_s2ccoin_fb"   , 600,  -15,  15);
  h_s2ccoin_fbac = new TH1F("h_s2ccoin_fbac" ,"h_s2ccoin_fbac" , 600,  -15,  15);
  h_s2ccoin_fbacz= new TH1F("h_s2ccoin_fbacz","h_s2ccoin_fbacz", 600,  -15,  15);
  h_s0ccoin_fb   = new TH1F("h_s0ccoin_fb"   ,"h_s0ccoin_fb"   , 600,-10,1000);
  h_Lz           = new TH1F("h_Lz"           ,"hLz"            , 200, -0.2, 0.2);
  h_Rz           = new TH1F("h_Rz"           ,"hRz"            , 200, -0.2, 0.2);
  set->SetTH1(h_s2coin_f1     ,"Coin time(S2 F1TDC)"                          ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s0coin_f1     ,"Coin time(S0 F1TDC)"                          ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s2ccoin_f1    ,"Coin time(S2 F1TDC w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s0ccoin_f1    ,"Coin time(S0 F1TDC w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s2ccoin_f1ac  ,"Coin time(S2 F1TDC w/ path cor. w/ AC cut)"   ,"coin time [ns]","counts", 2,3000,0);
  set->SetTH1(h_s2ccoin_f1acz ,"Coin time(S2 F1TDC w/ path cor. w/ AC&Z cut))","coin time [ns]","counts", 4,3000,0);
  set->SetTH1(h_s2ccoin_f1z   ,"Coin time(S2 F1TDC w/ path cor. w/ Z cut))"   ,"coin time [ns]","counts", 4,3000,0);
  set->SetTH1(h_s0ccoin_f1    ,"Coin time(S0 F1TDC w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s2coin_fb     ,"Coin time(S2 Fbus)"                          ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s0coin_fb     ,"Coin time(S0 Fbus w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s2ccoin_fb    ,"Coin time(S2 Fbus w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_s2ccoin_fbac  ,"Coin time(S2 Fbus w/ path cor. w/ AC cut)"   ,"coin time [ns]","counts", 2,3000,0);
  set->SetTH1(h_s2ccoin_fbacz ,"Coin time(S2 Fbus w/ path cor. w/ AC&Z cut))","coin time [ns]","counts", 4,3000,0);
  set->SetTH1(h_s0ccoin_fb    ,"Coin time(S0 Fbus w/ path cor.)"             ,"coin time [ns]","counts", 1,3000,0);
  set->SetTH1(h_Lz            ,"Target Z pos (L-HRS)"                        ,"z[m]"          ,"counts", 4,3000,0);
  set->SetTH1(h_Rz            ,"Target Z pos (R-HRS)"                        ,"z[m]"          ,"counts", 4,3000,0);


  h2_s2coin_f1_lpath     = new TH2F("h2_s2coin_f1_lpath"        ,"h2_s2coin_f1_lpath"        , 1000,   -15,   15,  100, 28.3,  29.8);
  h2_s0coin_f1_lpath     = new TH2F("h2_s0coin_f1_lpath"        ,"h2_s0coin_f1_lpath"        , 1000, -1000, 1000,  100,   25,  26.5);
  h2_s2coin_f1_rpath     = new TH2F("h2_s2coin_f1_rpath"        ,"h2_s2coin_f1_rpath"        , 1000,   -15,   15,  100, 28.5,  29.8);
  h2_s0coin_f1_rpath     = new TH2F("h2_s0coin_f1_rpath"        ,"h2_s0coin_f1_rpath"        , 1000, -1000, 1000,  100,   25,  26.5);
  h2_s2ccoin_f1_lpath    = new TH2F("h2_s2ccoin_f1_lpath"       ,"h2_s2ccoin_f1_lpath"       , 1000,   -15,   15,  100, 28.3,  29.8);
  h2_s0ccoin_f1_lpath    = new TH2F("h2_s0ccoin_f1_lpath"       ,"h2_s0ccoin_f1_lpath"       , 1000, -1000, 1000,  100,   25,  26.5);
  h2_s2ccoin_f1_rpath    = new TH2F("h2_s2ccoin_f1_rpath"       ,"h2_s2ccoin_f1_rpath"       , 1000,   -15,   15,  100, 28.3,  29.8);
  h2_s0ccoin_f1_rpath    = new TH2F("h2_s0ccoin_f1_rpath"       ,"h2_s0ccoin_f1_rpath"       , 1000, -1000, 1000,  100,   25,  26.5);
  h2_s2ccoin_f1_a1       = new TH2F("h2_s2ccoin_f1_a1"          ,"h2_s2ccoin_f1_ a1"         , 1000,   -15,   15,  100,-100., 4000.);
  h2_s2ccoin_f1_a2       = new TH2F("h2_s2ccoin_f1_a2"          ,"h2_s2ccoin_f1_a2"          , 1000,   -15,   15,  100,-100.,20000.);
  h2_s2ccoin_f1_a1_wa2cut= new TH2F("h2_s2ccoin_f1_a1_wa2cut"   ,"h2_s2ccoin_f1_a1_wa2cut"   , 1000,   -15,   15,  100,-100., 4000.);
  h2_s2ccoin_f1_a2_wa1cut= new TH2F("h2_s2ccoin_f1_a2_wa1cut"   ,"h2_s2ccoin_f1_a2_wa1cut"   , 1000,   -15,   15,  100,-100.,20000.);
  h2_s2coin_fb_lpath     = new TH2F("h2_s2coin_fb_lpath"        ,"h2_s2coin_fb_lpath"        ,  600,   -15,   15,  100, 28.3,  29.8);
  h2_s0coin_fb_lpath     = new TH2F("h2_s0coin_fb_lpath"        ,"h2_s0coin_fb_lpath"        , 1000, -1000, 1000,  100,   25,  26.5); 
  h2_s2coin_fb_rpath     = new TH2F("h2_s2coin_fb_rpath"        ,"h2_s2coin_fb_rpath"        ,  600,   -15,   15,  100, 28.5,  29.8);
  h2_s0coin_fb_rpath     = new TH2F("h2_s0coin_fb_rpath"        ,"h2_s0coin_fb_rpath"        , 1000, -1000, 1000,  100,   25,  26.5); 
  h2_s2ccoin_fb_lpath    = new TH2F("h2_s2ccoin_fb_lpath"       ,"h2_s2ccoin_fb_lpath"       ,  600,   -15,   15,  100, 28.3,  29.8);
  h2_s0ccoin_fb_lpath    = new TH2F("h2_s0ccoin_fb_lpath"       ,"h2_s0ccoin_fb_lpath"       , 1000, -1000, 1000,  100,   25,  26.5); 
  h2_s2ccoin_fb_rpath    = new TH2F("h2_s2ccoin_fb_rpath"       ,"h2_s2ccoin_fb_rpath"       ,  600,   -15,   15,  100, 28.3,  29.8);
  h2_s0ccoin_fb_rpath    = new TH2F("h2_s0ccoin_fb_rpath"       ,"h2_s0ccoin_fb_rpath"       , 1000, -1000, 1000,  100,   25,  26.5); 
  h2_s2ccoin_fb_a1       = new TH2F("h2_s2ccoin_fb_a1"          ,"h2_s2ccoin_fb_ a1"         ,  600,   -15,   15,  100,-100., 4000.);
  h2_s2ccoin_fb_a2       = new TH2F("h2_s2ccoin_fb_a2"          ,"h2_s2ccoin_fb_a2"          ,  600,   -15,   15,  100,-100.,20000.);
  h2_tof_beta            = new TH2F("h2_tof_beta"               ,"h2_tof_beta"               , 2000,  -100,  100,  200,   -2,     2);
  h2_msq_beta            = new TH2F("h2_msq_beta"               ,"h2_msq_beta"               , 2000,    -1,    1,  200,   -2,     2);
  h2_mom_beta            = new TH2F("h2_mom_beta"               ,"h2_mom_beta"               , 2000,   0.5,  3.5,  200, 1.72,  1.93);
  h2_path_beta           = new TH2F("h2_path_beta"              ,"h2_path_beta"              , 2000,    20,   30,  200,   -2,     2);
  h2_rf_path             = new TH2F("h2_rf_path"                ,"h2_rf_path"                , 2000,   -10,   10,  200,   25,  26.5);
  h2_LRz                 = new TH2F("h2_LRz"                    ,"h2_LRz"                    , 200,   -0.2,  0.2,  200, -0.2,   0.2);


  set->SetTH2(h2_s2coin_fb_lpath  ,"coin time(Fb) vs Path(LHRS)"        ,"cointime[ns]"    ,"Path(L-HRS)[m]");
  set->SetTH2(h2_s2coin_fb_rpath  ,"coin time(Fb) vs Path(RHRS)"        ,"cointime[ns]"    ,"Path(R-HRS)[m]");
  set->SetTH2(h2_s2ccoin_fb_lpath ,"coin time(Fb) vs Path(LHRS) w/ cor" ,"cointime[ns]"    ,"Path(L-HRS)[m]");
  set->SetTH2(h2_s2ccoin_fb_rpath ,"coin time(Fb) vs Path(RHRS) w/ cor" ,"cointime[ns]"    ,"Path(R-HRS)[m]");
  set->SetTH2(h2_s2ccoin_fb_a1    ,"coin time(Fb) vs A1 ADC w/ cor"     ,"cointime[ns]"    ,"A1 ADCp[arb]");
  set->SetTH2(h2_s2ccoin_fb_a2    ,"coin time(Fb) vs A1 ADC w/ cor"     ,"cointime[ns]"    ,"A2 ADCp[arb]");

  set->SetTH2(h2_s2coin_f1_lpath  ,"coin time(F1) vs Path(LHRS)"        ,"cointime[ns]"    ,"Path(L-HRS)[m]");
  set->SetTH2(h2_s2coin_f1_rpath  ,"coin time(F1) vs Path(RHRS)"        ,"cointime[ns]"    ,"Path(R-HRS)[m]");
  set->SetTH2(h2_s2ccoin_f1_lpath ,"coin time(F1) vs Path(LHRS) w/ cor" ,"cointime[ns]"    ,"Path(L-HRS)[m]");
  set->SetTH2(h2_s2ccoin_f1_rpath ,"coin time(F1) vs Path(RHRS) w/ cor" ,"cointime[ns]"    ,"Path(R-HRS)[m]");
  set->SetTH2(h2_s2ccoin_f1_a1    ,"coin time(F1) vs A1 ADC w/ cor"     ,"cointime[ns]"    ,"A1 ADCp[arb]");
  set->SetTH2(h2_s2ccoin_f1_a2    ,"coin time(F1) vs A2 ADC w/ cor"     ,"cointime[ns]"    ,"A2 ADCp[arb]");
  set->SetTH2(h2_s2ccoin_f1_a1_wa2cut    ,"ccoin time(F1) vs A1 ADC w/ A2cut"     ,"cointime[ns]"    ,"A1 ADCp[arb]");
  set->SetTH2(h2_s2ccoin_f1_a2_wa1cut    ,"ccoin time(F1) vs A2 ADC w/ A1cut"     ,"cointime[ns]"    ,"A2 ADCp[arb]");
  
  set->SetTH2(h2_tof_beta ,"ToF vs beta" ,"ToF"    ,"#beta");
  set->SetTH2(h2_msq_beta ,"Msq vs beta" ,"m^{2}"  ,"#beta");
  set->SetTH2(h2_mom_beta ,"Mom vs beta" ,"1/#beta","mom"  );
  set->SetTH2(h2_path_beta,"Path vs beta","pathl"  ,"#beta");
  set->SetTH2(h2_rf_path  ,"RF-S2 vs Pathl","RF-S2 [ns]"  ,"path[m]");
  set->SetTH2(h2_LRz      ,"L-HRS z vs R-HRS z"     ,"Lz[m]"    ,"Rz[m]");

  for(int i=0;i<RS2;i++){
    h_s2r_tof[i]  = new TH1F(Form("h_s2r_tof",i) , Form("h_s2r_tof",i+1)      ,2000,-100,100);
    h_s2r_msq[i]  = new TH1F(Form("h_s2r_msq",i) , Form("h_s2r_msq",i+1)      , 200,  -1,  1);
    h_s2r_beta[i] = new TH1F(Form("h_s2r_beta",i), Form("h_s2r_beta",i+1)     , 100,   0,  2);
    set->SetTH1(h_s2r_tof[i]  ,Form("ToF(S2R%d - S0R)",i),"time[ns]","counts");
    set->SetTH1(h_s2r_msq[i]  ,Form("Mass^2F(S2R%d)",i)  ,"mass[]","counts");
    set->SetTH1(h_s2r_beta[i] ,Form("#beta(S2R%d)",i)    ,"#beta","counts");

    h2_s2r_tof_beta[i]    = new TH2F(Form("h2_s2r_tof_beta%d", i+1)  , Form("h2_s2r_tof_beta%d", i) , 2000, -100,  100,  200,   -2,   2);
    h2_s2r_msq_beta[i]    = new TH2F(Form("h2_s2r_msq_beta%d", i+1)  , Form("h2_s2r_msq_beta%d", i) , 2000,   -1,    1,  200,   -2,   2);
    h2_s2r_mom_beta[i]    = new TH2F(Form("h2_s2r_mom_beta%d", i+1)  , Form("h2_s2r_mom_beta%d", i) , 2000,  -.5,  2.5,  200,  1.6,   2);
    h2_s2r_path_beta[i]   = new TH2F(Form("h2_s2r_path_beta%d", i+1) , Form("h2_s2r_path_beta%d", i), 2000,   20,   30,  200,   -2,   2);
    set->SetTH2(h2_s2r_mom_beta[i]  ,Form("p vs 1/#beta(S2R%d - S0R)",i),"1/#beta","mom[GeV/c]");

  }

  //ac1 study (pion selection by cointime cut)
  h2_a1_asump_x_pi = new TH2F("h2_a1_asump_x_pi","h2_a1_asump_x_pi", 400,   -1000,  10000,  400, -1. ,   1. );
  h2_a1_asump_y_pi = new TH2F("h2_a1_asump_y_pi","h2_a1_asump_y_pi", 400,   -1000,  10000,  400, -0.2,   0.2);
  set->SetTH2(h2_a1_asump_x_pi,"a1 adc sum vs x[m]","a1 adc sum","x[m]");
  set->SetTH2(h2_a1_asump_y_pi,"a1 adc sum vs y[m]","a1 adc sum","y[m]");

  for(int i=0;i<nhist_a1_asump_pi;i++){
    h_a1_asump_pi[i] = new TH1D(Form("h_a1_asump_pi%d",i+1),Form("h_a1_asump_pi%d",i+1), 400,   -1000,  10000);
    set->SetTH1(h_a1_asump_pi[i],Form("h_a1_asump_pi%d",i+1),"a1 adc sum","counts");
  }

}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::loop(){

  bool Rtr_flag, Ltr_flag;
  bool RZ_flag, LZ_flag, Z_flag;
  
  int eventcount=0;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%10000==0)cout<<n <<" / "<<ENum<<endl;
    tree->GetEntry(n);
    convertF1TDCR(param);
    convertF1TDCL(param);

    double Rm2[MAX],RF_S2[MAX],beta_k[MAX], Vk[MAX];
    int ritrack = 0, litrack = 0;

    if(DR_T5>0.){//coin trig.

#ifdef ALL_TRACK
   for(ritrack=0;ritrack<R_tr_n;ritrack++){
     for(litrack=0;litrack<L_tr_n;litrack++){
#endif

      Rtr_flag=Ltr_flag=RZ_flag=LZ_flag=false;
      ///R-HRS////
      Rm2[ritrack] = ( 1. / (R_tr_beta[ritrack] * R_tr_beta[ritrack]) -1. ) * R_tr_p[ritrack] * R_tr_p[ritrack];
      beta_k[ritrack] = R_tr_p[ritrack]/sqrt(R_tr_p[ritrack]*R_tr_p[ritrack] + MK*MK);
      Vk[ritrack] = beta_k[ritrack]*LightVelocity;


      h_Rz -> Fill(R_tr_vz[ritrack]);
      h2_msq_beta   ->Fill(Rm2[ritrack]         , R_tr_beta[ritrack]);
      h2_mom_beta   ->Fill(1./R_tr_beta[ritrack], R_tr_p[ritrack]   );
      h2_path_beta  ->Fill(R_tr_pathl[ritrack]  , R_tr_beta[ritrack]);
      if(R_s2_trpad[ritrack]==8)h2_rf_path    ->Fill(RF_S2[ritrack]       , R_tr_pathl[ritrack]   );
      for(int i=0;i<RS2;i++){
        if(R_s2_trpad[ritrack]==i  && R_s0_lt[0]>0&& R_s0_rt[0]>0){
          h2_tof_beta   ->Fill(s2ns*(R_s2_time[i] - R_s0_time[0]) , R_tr_beta[ritrack]);
          h_s2r_tof[i]  ->Fill(s2ns*(R_s2_time[i] - R_s0_time[0]));
          h_s2r_msq[i]  ->Fill(Rm2[ritrack]);
          h_s2r_beta[i] ->Fill(R_tr_beta[ritrack]);

          h2_s2r_tof_beta[i]  ->Fill(s2ns*(R_s2_time[i] - R_s0_time[0]) , R_tr_beta[ritrack]);
          h2_s2r_msq_beta[i]  ->Fill(Rm2[ritrack]                       , R_tr_beta[ritrack]);
          h2_s2r_mom_beta[i]  ->Fill(1./R_tr_beta[ritrack]              , R_tr_p[ritrack]   );
          h2_s2r_path_beta[i] ->Fill(R_tr_pathl[ritrack]                , R_tr_beta[ritrack]);
        }
      }//each S2 paddle

      ///L-HRS////

      h_Lz -> Fill(L_tr_vz[litrack]);
      h2_LRz -> Fill(L_tr_vz[litrack],R_tr_vz[ritrack]);

    if(R_tr_pathl[ritrack]+R_s2_trpath[ritrack]>28.7 && R_tr_pathl[ritrack]+R_s2_trpath[ritrack]<29.4)Rtr_flag=true;
    if(L_tr_pathl[litrack]+L_s2_trpath[litrack]>28.6 && L_tr_pathl[litrack]+L_s2_trpath[litrack]<29.2)Ltr_flag=true;
    if(R_tr_vz[ritrack]>-0.1 && R_tr_vz[ritrack]< 0.1)RZ_flag=true;
    if(L_tr_vz[litrack]>-0.1 && L_tr_vz[litrack]< 0.1)LZ_flag=true;

    double f1s2coin  =  -1.*(LS2_F1time[(int)L_s2_trpad[litrack]] - RS2_F1time[(int)R_s2_trpad[ritrack]] );
    double f1s2ccoin = -1.*((LS2_F1time[(int)L_s2_trpad[litrack]] - L_tr_pathl[litrack]/LightVelocity) - (RS2_F1time[(int)R_s2_trpad[ritrack]] - R_tr_pathl[ritrack]/Vk[ritrack]) );
    f1s2coin  -= f1offset;
    f1s2ccoin -= f1coffset;

    double fbs2coin  = s2ns*(L_s2_time[(int)L_s2_trpad[litrack]] - R_s2_time[(int)R_s2_trpad[ritrack]]);
    double fbs2ccoin = (s2ns*L_s2_time[(int)L_s2_trpad[litrack]] + (L_tr_pathl[litrack]+L_s2_trpath[litrack])/LightVelocity ) - (s2ns*R_s2_time[(int)R_s2_trpad[ritrack]] + (R_tr_pathl[ritrack]+R_s2_trpath[ritrack])/Vk[ritrack]);

    fbs2coin  -= coffset;
    fbs2ccoin -= ccoffset;
    //double fbs2ccoin = (s2ns*L_s2_time[(int)L_s2_trpad[litrack]] + L_tr_pathl[litrack]/LightVelocity) - (s2ns*R_s2_time[(int)R_s2_trpad[ritrack]] - R_tr_pathl[ritrack]/Vk[ritrack]);
      ///cointime L-R///
    if(R_s2_rt[(int)R_s2_trpad[ritrack]]>1. && R_s2_lt[(int)R_s2_trpad[ritrack]]>1. && L_s2_rt[(int)L_s2_trpad[litrack]]>1.&& L_s2_lt[(int)L_s2_trpad[litrack]]>1.){
      if(Rtr_flag && Ltr_flag){
        h_s2coin_fb          ->Fill(fbs2coin);
        h2_s2coin_fb_lpath   ->Fill(fbs2coin ,  L_tr_pathl[litrack]+L_s2_trpath[litrack]);
        h2_s2coin_fb_rpath   ->Fill(fbs2coin ,  R_tr_pathl[ritrack]+R_s2_trpath[ritrack]);
        h_s2ccoin_fb         ->Fill(fbs2ccoin);
        h2_s2ccoin_fb_lpath  ->Fill(fbs2ccoin,  L_tr_pathl[litrack]+L_s2_trpath[litrack]);
        h2_s2ccoin_fb_rpath  ->Fill(fbs2ccoin,  R_tr_pathl[ritrack]+R_s2_trpath[ritrack]);
        h2_s2ccoin_fb_a1     ->Fill(fbs2ccoin,  R_a1_asum_p);
        h2_s2ccoin_fb_a2     ->Fill(fbs2ccoin,  R_a2_asum_p);
        if(R_a1_asum_p<150. && R_a2_asum_p>1400.)h_s2ccoin_fbac    ->Fill(fbs2ccoin);
        if(R_a1_asum_p<150. && R_a2_asum_p>1400. && RZ_flag && LZ_flag)h_s2ccoin_fbacz   ->Fill(fbs2ccoin);
        //cout<<  R_a2_asum_p <<endl;
      }
    }
    if(RS2T_F1TDC[(int)R_s2_trpad[ritrack]]>1. && RS2B_F1TDC[(int)R_s2_trpad[ritrack]]>1. && LS2T_F1TDC[(int)L_s2_trpad[litrack]]>1.&& LS2B_F1TDC[(int)L_s2_trpad[litrack]]>1.){
      if(Rtr_flag && Ltr_flag){
        h_s2coin_f1  ->Fill(f1s2coin);
        h_s2ccoin_f1 ->Fill(f1s2ccoin);
        h2_s2coin_f1_lpath   ->Fill(f1s2coin ,  L_tr_pathl[litrack]+L_s2_trpath[litrack]);
        h2_s2coin_f1_rpath   ->Fill(f1s2coin ,  R_tr_pathl[ritrack]+R_s2_trpath[ritrack]);
        h_s2ccoin_f1         ->Fill(f1s2ccoin);
        h2_s2ccoin_f1_lpath  ->Fill(f1s2ccoin,  L_tr_pathl[litrack]+L_s2_trpath[litrack]);
        h2_s2ccoin_f1_rpath  ->Fill(f1s2ccoin,  R_tr_pathl[ritrack]+R_s2_trpath[ritrack]);
        h2_s2ccoin_f1_a1     ->Fill(f1s2ccoin,  R_a1_asum_p);
        h2_s2ccoin_f1_a2     ->Fill(f1s2ccoin,  R_a2_asum_p);
        if(R_a1_asum_p<150. && R_a2_asum_p>1400.)h_s2ccoin_f1ac    ->Fill(f1s2ccoin);
        if(R_a1_asum_p<150. && R_a2_asum_p>1400. && RZ_flag && LZ_flag)h_s2ccoin_f1acz   ->Fill(f1s2ccoin);
        if(R_a2_asum_p>1400. && RZ_flag && LZ_flag)h2_s2ccoin_f1_a1_wa2cut   ->Fill(f1s2ccoin, R_a1_asum_p);
        if(R_a1_asum_p<150.  && RZ_flag && LZ_flag)h2_s2ccoin_f1_a2_wa1cut   ->Fill(f1s2ccoin, R_a2_asum_p);
      }
    }

    if(R_s0_rt[0]>1. && R_s0_lt[0]>1. && L_s0_rt[0]>1.&& L_s0_lt[0]>1. ){

      if(LF1Ref[0]>1. && LF1Ref[1]>1. && RF1Ref[0]>1. && RF1Ref[1]>1.){
        //cout<<"LF1Ref : "<<LF1Ref[0]<< " "<<LF1Ref[1]<<", RF1Ref : "<<RF1Ref[0]<<" "<<RF1Ref[1]<<endl;
        h_s0coin_f1  ->Fill(LS0_F1time[0] - RS0_F1time[0]);
        h_s0ccoin_f1 ->Fill(LS0_F1time[0] - RS0_F1time[0]);
      }

      h_s0coin_fb  ->Fill(s2ns*(L_s2_time[0] - R_s2_time[0]));
      h_s0ccoin_fb ->Fill(LS0_F1time[0] - RS0_F1time[0]);

    }

    //ac1 study (pion selection by cointime cut)
    if(RS2T_F1TDC[(int)R_s2_trpad[ritrack]]>1. && RS2B_F1TDC[(int)R_s2_trpad[ritrack]]>1. && LS2T_F1TDC[(int)L_s2_trpad[litrack]]>1.&& LS2B_F1TDC[(int)L_s2_trpad[litrack]]>1.){
      if(Rtr_flag && Ltr_flag){
        if(RZ_flag && LZ_flag)h_s2ccoin_f1z   ->Fill(f1s2ccoin);
        if(RZ_flag && LZ_flag && f1s2ccoin>-1. && f1s2ccoin<1.){
	  h2_a1_asump_x_pi ->Fill(R_a1_asum_p, R_a1_trx[0]);
	  h2_a1_asump_y_pi ->Fill(R_a1_asum_p, R_a1_try[0]);
	  // eventcount++;
	}
      }
    }

      
#ifdef ALL_TRACK
   }}
#endif
    }//coin_trig
  }
  // cout<<eventcount<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::fit(){
  h_frame[0] = new TH2F("h_frame0","h_frame0",10, -0.5, 15.5,10,0.8,1.2);
  h_frame[1] = new TH2F("h_frame1","h_frame1",10, -0.5, 15.5,10,0.0,0.3);
  set->SetTH2(h_frame[0] , "#beta peak pos each S2","S2 paddle","#beta peak");
  set->SetTH2(h_frame[1] , "#beta width each S2"   ,"S2 paddle","#beta width");
  tg_beta_pos = new TGraphErrors(); 
  tg_beta_wid = new TGraphErrors();
  set->SetGrErr(tg_beta_pos, "#beta peak pos each S2","S2 paddle","#beta peak",1,4,23);
  set->SetGrErr(tg_beta_wid, "#beta width each S2"   ,"S2 paddle","#beta width",1,4,24);

  for(int i=0;i<16;i++){
    ga_beta[i] = new TF1(Form("ga_beta%d",i+1),"gaus",-2,2);
    set->SetTF1(ga_beta[i],2,1,1);
    double min=0.5,max=1.5;
    set->FitGaus(h_s2r_beta[i],min,max,1.5,5);
    h_s2r_beta[i]->Fit(ga_beta[i],"QR","",min,max);
    beta_pos[i]  = ga_beta[i]->GetParameter(1);
    beta_wid[i]  = ga_beta[i]->GetParameter(2);
    ebeta_pos[i] = ga_beta[i]->GetParError(1);
    ebeta_wid[i] = ga_beta[i]->GetParError(2);;

    tg_beta_pos ->SetPoint(i,i,beta_pos[i]); 
    tg_beta_wid ->SetPoint(i,i,beta_wid[i]);
    tg_beta_pos ->SetPointError(i,0,ebeta_pos[i]); 
    tg_beta_wid ->SetPointError(i,0,ebeta_wid[i]);
    
  }

  
  ga_coin = new TF1("ga_coin","gaus",-100,500);
  set->SetTF1(ga_coin,2,1,1);
  double min=-40.0,max= 50.0;
  set->FitGaus(h_s2ccoin_f1,min,max,1.5,5);
  h_s2ccoin_f1->Fit(ga_coin,"QR","",min,max);

  cout<<"cointime peak(pion)= "<<ga_coin->GetParameter(1)<<", width= "<<ga_coin->GetParameter(2)<<" (ns)"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::projection(){

  int bin_threshold;
  int threshold_a1_pi = 150;
  double integral_all,integral_threshold;
  double a1_eff[nhist_a1_asump_pi];
  double bin_mean;
  double x_mean[nhist_a1_asump_pi];

  for(int i=0;i<nhist_a1_asump_pi;i++){
    h_a1_asump_pi[i] = (TH1D*)h2_a1_asump_x_pi->ProjectionX(Form("h_a1_asump_pi%d",i+1),i*nbin_a1_asump_pi+1,(i+1)*nbin_a1_asump_pi,"");
    bin_threshold = h_a1_asump_pi[i]->GetXaxis()->FindBin(threshold_a1_pi);
    integral_all = h_a1_asump_pi[i]->Integral(1,400,"");
    integral_threshold = h_a1_asump_pi[i]->Integral(bin_threshold,400,"");
    a1_eff[i] = integral_threshold / integral_all;
    bin_mean = i*nbin_a1_asump_pi + 0.5*nbin_a1_asump_pi;
    x_mean[i] = h2_a1_asump_x_pi->GetYaxis()->GetBinCenter(bin_mean);
  }
  g_a1_eff_pi = new TGraph(nhist_a1_asump_pi,x_mean,a1_eff);
  set->SetGr(g_a1_eff_pi,"A1_Efficiency_x-pos","x-position [m]","Efficiency",2,2,20,0);
  g_a1_eff_pi->Write();
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::draw(){

  //F1
  c[0]->Clear();c[0]->Divide(4,2);
  c[0]->cd(1);h_s2coin_f1  ->Draw();
  c[0]->cd(2);h_s2ccoin_f1 ->Draw();
  c[0]->cd(3);gPad->SetLogz(1);h2_s2ccoin_f1_a1  ->Draw("colz");
  c[0]->cd(4);gPad->SetLogz(1);h2_s2ccoin_f1_a2  ->Draw("colz");
  c[0]->cd(5);h2_s2coin_f1_lpath   ->Draw("colz");
  c[0]->cd(6);h2_s2coin_f1_rpath   ->Draw("colz");
  c[0]->cd(7);h2_s2ccoin_f1_lpath  ->Draw("colz");
  c[0]->cd(8);h2_s2ccoin_f1_rpath  ->Draw("colz");
  //Fbus
  c[1]->Clear();c[1]->Divide(4,2);
  c[1]->cd(1);h_s2coin_fb  ->Draw();
  c[1]->cd(2);h_s2ccoin_fb ->Draw();
  c[1]->cd(3);gPad->SetLogz(1);h2_s2ccoin_fb_a1  ->Draw("colz");
  c[1]->cd(4);gPad->SetLogz(1);h2_s2ccoin_fb_a2  ->Draw("colz");
  c[1]->cd(5);h2_s2coin_fb_lpath   ->Draw("colz");
  c[1]->cd(6);h2_s2coin_fb_rpath   ->Draw("colz");
  c[1]->cd(7);h2_s2ccoin_fb_lpath  ->Draw("colz");
  c[1]->cd(8);h2_s2ccoin_fb_rpath  ->Draw("colz");

  c[2]->Clear();c[2]->Divide(1,1);
  c[2]->cd(1);h_s2ccoin_f1ac ->Draw();
  h_s2ccoin_f1ac_rebin = (TH1F*)h_s2ccoin_f1ac ->Clone();
  h_s2ccoin_f1ac_rebin ->Rebin(4);
  //c[2]->cd(1);h_s2ccoin_fbac ->Draw();
  //h_s2ccoin_fbac_rebin = (TH1F*)h_s2ccoin_fbac ->Clone();
  //h_s2ccoin_fbac_rebin ->Rebin(4);

  c[3]->Clear();c[3]->Divide(2,2);
  c[3]->cd(1);h_s2ccoin_f1  ->Draw();h_s2ccoin_f1ac ->Draw("same");
  c[3]->cd(2);gPad->SetLogy(1);h_s2ccoin_f1  ->Draw();h_s2ccoin_f1ac ->Draw("same");
  c[3]->cd(3);gPad->SetLogy(1);h_s2ccoin_f1  ->Draw();h_s2ccoin_f1ac_rebin ->Draw("same");
  c[3]->cd(4);gPad->SetLogy(0);h_s2ccoin_f1ac_rebin ->Draw("");
  //c[3]->cd(1);h_s2ccoin_fb  ->Draw();h_s2ccoin_fbac ->Draw("same");
  //c[3]->cd(2);gPad->SetLogy(1);h_s2ccoin_fb  ->Draw();h_s2ccoin_fbac ->Draw("same");
  //c[3]->cd(3);gPad->SetLogy(1);h_s2ccoin_fb  ->Draw();h_s2ccoin_fbac_rebin ->Draw("same");
  //c[3]->cd(4);gPad->SetLogy(0);h_s2ccoin_fbac_rebin ->Draw("");



  c[4]->Clear();c[4]->Divide(3,3);
  c[4]->cd(1);gPad->SetLogz(1);h2_msq_beta   ->Draw("colz");
  c[4]->cd(2);gPad->SetLogz(1);h2_mom_beta   ->Draw("colz");
  c[4]->cd(3);h_Lz   -> Draw("");
  c[4]->cd(4);h_Rz   -> Draw("");
  c[4]->cd(5);h2_LRz -> Draw("colz");
  c[4]->cd(6);h_s2ccoin_fbacz ->Draw();
  c[4]->cd(7);gPad->SetLogz(1);h2_s2ccoin_f1_a1_wa2cut  ->Draw("colz");
  c[4]->cd(8);gPad->SetLogz(1);h2_s2ccoin_f1_a2_wa1cut  ->Draw("colz");

  c[5]->Clear();c[5]->Divide(1,1);
  c[5]->cd(1);h_s2ccoin_f1->Draw();
  //c[4]->cd(3);gPad->SetLogz(1);h_frame[0]->Draw();tg_beta_pos ->Draw("sameP");//h2_path_beta  ->Draw("colz");
  //c[4]->cd(4);gPad->SetLogz(1);h_frame[1]->Draw();tg_beta_wid ->Draw("sameP");//h2_tof_beta   ->Draw("colz");

  //c[3]->Clear();c[3]->Divide(4,4);
  //for(int i=0;i<16;i++){
  //  c[3]->cd(i+1);gPad->SetLogz(1);h2_s2r_tof_beta[i]  ->Draw("colz");
  //}
  //c[5]->Clear();//c[5]->Divide(2,2);
  //c[5]->cd(1);gPad->SetLogz(1);h2_rf_path->Draw("colz");
}
////////////////////////////////////////////////////////////////////////////
void ana_cointime::savecanvas(string ofname){
  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[NCanvas-1]->Print(Form("%s]",ofname.c_str()) );
  cout<<ofname<<" saved"<<endl;
  ofp->Write();
  ofp->Close();
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/cosmic1020.root";
  string runlistname = "runlist/test.txt";
  string ofname = "tmp.root";
  string ofname_pdf = "tmp.pdf";
  string paramname = "twlk_param/default.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool single_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:s:w:n:bop:"))!=-1){
    switch(ch){
    case 'f':
      runlistname = optarg;
      cout<<"input runlist filename : "<<runlistname<<endl;
      break;
    case 's':
      single_flag = true;
      ifname = optarg;
      cout<<"input root filename : "<<ifname<<endl;
      cout<<"signle file analysis mode"<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      draw_flag = false;
      paramname = optarg;
      cout<<"input param name : "<<paramname<<endl;
      break;
    case 'h':
      cout<<"-f : input runlist file"<<endl;
      cout<<"-s : input root file name(single run analysis)"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : input param file"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  TApplication *theApp = new TApplication("App", &argc, argv);
  ana_cointime *coin = new ana_cointime();


  ofname_pdf = ofname;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append(".pdf");

  coin->SetMaxEvent(MaxNum);
  coin->SetInputParam(paramname);
  coin->SetOutputRoot(ofname);
  if(single_flag)coin->SetRoot(ifname);
  else coin-> SetRunList(runlistname);
  coin->makehist();
  coin->loop();
  coin->fit();
  coin->projection();
  coin->draw();
  if(output_flag)coin->savecanvas(ofname_pdf);
  delete coin;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

