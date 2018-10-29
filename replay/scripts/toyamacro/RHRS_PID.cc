////////////////////////////
//PID plots HRSR          //
// by Y. Toyama Oct. 2018 //
////////////////////////////

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

const int NCanvas = 5;//num of canvas

class RHRS_PID : public Tree
{
 public:
         RHRS_PID();
        ~RHRS_PID();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void SetInputParam(string ifname);
  Setting *set;
  ParamMan *param;

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;

    //General
    TH1F *h_tof, *h_beta, *h_msq, *h_mom, *h_path, *h_t_at_targ;
    TH2F *h2_tof_beta, *h2_beta, *h2_msq_beta, *h2_mom_beta, *h2_path_beta;

    //Each S2 paddle
    TH1F *h_s2r_tof[RS2], *h_s2r_beta[RS2], *h_s2r_msq[RS2];
    TH2F *h2_s2r_tof_beta[RS2], *h2_s2r_msq_beta[RS2], *h2_s2r_mom_beta[RS2], *h2_s2r_path_beta[RS2];

    int run_num;
    TCanvas *c[NCanvas];

    TF1 *ga_beta[16];
    TGraphErrors *tg_beta_pos, *tg_beta_wid;
    TH2F *h_frame[2];

    double beta_pos[16],beta_wid[16];
    double ebeta_pos[16],ebeta_wid[16];
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RHRS_PID::RHRS_PID()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(0);

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
RHRS_PID::~RHRS_PID(){
}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::SetInputParam(string ifname){
  param = new ParamMan(ifname.c_str());
  if(param -> SetVal())cout<<"F1TDC parameter setted"<<endl; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::makehist(){
  h2_tof_beta    = new TH2F("h2_tof_beta"  ,"h2_tof_beta"  , 2000, -100,  100,  200,   -2,   2);
  h2_msq_beta    = new TH2F("h2_msq_beta"  ,"h2_msq_beta"  , 2000,   -1,    1,  200,   -2,   2);
  h2_mom_beta    = new TH2F("h2_mom_beta"  ,"h2_mom_beta"  , 2000, 0.5,  3.5,  200,   -2,   2);
  h2_path_beta   = new TH2F("h2_path_beta" ,"h2_path_beta" , 2000,  20,   30,  200,   -2,   2);
  set->SetTH2(h2_tof_beta ,"ToF vs beta" ,"ToF"  ,"#beta");
  set->SetTH2(h2_msq_beta ,"Msq vs beta" ,"m^{2}","#beta");
  set->SetTH2(h2_mom_beta ,"Mom vs beta" ,"mom"  ,"#beta");
  set->SetTH2(h2_path_beta,"Path vs beta","pathl","#beta");

  for(int i=0;i<RS2;i++){
    h_s2r_tof[i]  = new TH1F(Form("h_s2r_tof",i) , Form("h_s2r_tof",i+1)      ,2000,-100,100);
    h_s2r_msq[i]  = new TH1F(Form("h_s2r_msq",i) , Form("h_s2r_msq",i+1)      , 200,  -1,  1);
    h_s2r_beta[i] = new TH1F(Form("h_s2r_beta",i), Form("h_s2r_beta",i+1)     , 100,   0,  2);
    set->SetTH1(h_s2r_tof[i]  ,Form("ToF(S2R%d - S0R)",i),"time[ns]","counts");
    set->SetTH1(h_s2r_msq[i]  ,Form("Mass^2F(S2R%d)",i)  ,"mass[]","counts");
    set->SetTH1(h_s2r_beta[i] ,Form("#beta(S2R%d)",i)    ,"#beta","counts");

    h2_s2r_tof_beta[i]    = new TH2F(Form("h2_s2r_tof_beta%d", i+1)  , Form("h2_s2r_tof_beta%d", i) , 2000, -100,  100,  200,   -2,   2);
    h2_s2r_msq_beta[i]    = new TH2F(Form("h2_s2r_msq_beta%d", i+1)  , Form("h2_s2r_msq_beta%d", i) , 2000,   -1,    1,  200,   -2,   2);
    h2_s2r_mom_beta[i]    = new TH2F(Form("h2_s2r_mom_beta%d", i+1)  , Form("h2_s2r_mom_beta%d", i) , 2000,  0.5,  3.5,  200,   -2,   2);
    h2_s2r_path_beta[i]   = new TH2F(Form("h2_s2r_path_beta%d", i+1) , Form("h2_s2r_path_beta%d", i), 2000,   20,   30,  200,   -2,   2);

  }

}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%1000==0)cout<<n <<" / "<<ENum<<endl;
    tree->GetEntry(n);
    convertF1TDCR(param);

    double Rm2[MAX];
    for(int itrack=0;itrack<R_tr_n;itrack++){//each track
      if(DR_T4>0.){
        Rm2[itrack] = ( 1. / (R_tr_beta[itrack] * R_tr_beta[itrack]) -1. ) * R_tr_p[itrack] * R_tr_p[itrack];
        h2_msq_beta   ->Fill(Rm2[itrack]        , R_tr_beta[itrack]);
        h2_mom_beta   ->Fill(R_tr_p[itrack]*10. , R_tr_beta[itrack]);
        h2_path_beta  ->Fill(R_tr_pathl[itrack] , R_tr_beta[itrack]);
        for(int i=0;i<RS2;i++){
          if(R_s2_trpad[itrack]==i  && R_s0_lt[0]>0&& R_s0_rt[0]>0){
            h2_tof_beta   ->Fill(1.e+9*(R_s2_time[i] - R_s0_time[0]) , R_tr_beta[itrack]);
            h_s2r_tof[i]  ->Fill(1.e+9*(R_s2_time[i] - R_s0_time[0]));
            h_s2r_msq[i]  ->Fill(Rm2[itrack]);
            h_s2r_beta[i] ->Fill(R_tr_beta[itrack]);

            h2_s2r_tof_beta[i]  ->Fill(s2ns*(R_s2_time[i] - R_s0_time[0]), R_tr_beta[itrack]);
            h2_s2r_msq_beta[i]  ->Fill(Rm2[itrack]       , R_tr_beta[itrack]);
            h2_s2r_mom_beta[i]  ->Fill(R_tr_p[itrack]*10., R_tr_beta[itrack]);
            h2_s2r_path_beta[i] ->Fill(R_tr_pathl[itrack], R_tr_beta[itrack]);
          }
        }//each S2 paddle
      }//S0&&S2
    }//each track
  }

}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::fit(){
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
}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::draw(){

  c[0]->Clear();c[0]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[0]->cd(i+1);gPad->SetLogy(1);h_s2r_tof[i]->Draw();
  }

  c[1]->Clear();c[1]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[1]->cd(i+1);gPad->SetLogy(1);h_s2r_beta[i]->Draw();
  }

  c[2]->Clear();c[2]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[2]->cd(i+1);gPad->SetLogy(1);h_s2r_msq[i]->Draw();
  }
  c[3]->Clear();c[3]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[3]->cd(i+1);gPad->SetLogz(1);h2_s2r_tof_beta[i]  ->Draw("colz");
  }

  c[4]->Clear();c[4]->Divide(2,2);
  c[4]->cd(1);gPad->SetLogz(1);h2_msq_beta   ->Draw("colz");
  c[4]->cd(2);gPad->SetLogz(1);h2_mom_beta   ->Draw("colz");
  c[4]->cd(3);gPad->SetLogz(1);h_frame[0]->Draw();tg_beta_pos ->Draw("sameP");//h2_path_beta  ->Draw("colz");
  c[4]->cd(4);gPad->SetLogz(1);h_frame[1]->Draw();tg_beta_wid ->Draw("sameP");//h2_tof_beta   ->Draw("colz");

}
////////////////////////////////////////////////////////////////////////////
void RHRS_PID::savecanvas(string ofname){
  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[NCanvas-1]->Print(Form("%s]",ofname.c_str()) );
  cout<<ofname<<" saved"<<endl;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/cosmic1020.root";
  string ofname = "toyamacro/pdf/RHRS_PID1020.pdf";
  string paramname = "twlk_param/default.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bcop:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
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
    case 'c':
      coin_flag = true;
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
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
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
  RHRS_PID *pid = new RHRS_PID();

  pid->SetMaxEvent(MaxNum);
  pid->SetInputParam(paramname);
  pid->SetRoot(ifname);
  pid->makehist();
  pid->loop();
  pid->fit();
  pid->draw();
  if(output_flag)pid->savecanvas(ofname);
  delete pid;

  //gSystem->Exit(1);
  theApp->Run();
  return 0;
}

