//////////////////////////////////////
//t0 calibration of S2(F1TDC & Fbus)//
// by Y. Toyama Oct. 2018           //
//////////////////////////////////////

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

#define Calibration

const int NCanvas = 9;//num of canvas

class s2_t0_calib : public Tree
{
 public:
         s2_t0_calib();
        ~s2_t0_calib();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void SetInputParam(string ifname);
  void SetLR(int lr){LR=lr;}

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    bool anaL_oneevent();
    bool anaR_oneevent();

    TH1F *h_s2s0_beta[16], *h_s2s0_beta_FB[16];
    TH1F *h_s2s0_tof[16], *h_s2s0_tdiff[16], *h_s2s0_tof_FB[16], *h_s2s0_tdiff_FB[16];
    TH2F *h2_s2s0_f1beta_fbbeta[16], *h2_s2s0_f1tof_fbtof[16];
    TF1 *ga_tdiffF1[16], *ga_tdiffFB[16];
    TH2F *h_frame[4];
    int LR;//L = 0, R = 1
    int run_num;
    Setting *set;
    ParamMan *param;
    TCanvas *c[NCanvas];
    TGraphErrors *tg_tdiffFB_pos, *tg_tdiffFB_wid;
    TGraphErrors *tg_tdiffF1_pos, *tg_tdiffF1_wid;

    double tdiffFB_pos[16], tdiffFB_wid[16], etdiffFB_pos[16], etdiffFB_wid[16];
    double tdiffF1_pos[16], tdiffF1_wid[16], etdiffF1_pos[16], etdiffF1_wid[16];



    int tr_n;//num. of track
    double beta[MAX],betaF1[MAX],s2_trpad[MAX],paths2s0[MAX];
    double S2T_F1TDC[16],S0T_F1TDC[1];
    double S2B_F1TDC[16],S0B_F1TDC[1];
    double S2_F1time[16],S0_F1time[1];
    double S2_lt[16],S2_rt[16],S0_lt[1],S0_rt[1];
    double S2_time[16],S0_time[1];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s2_t0_calib::s2_t0_calib()
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
s2_t0_calib::~s2_t0_calib(){
}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  if(LR==0)     readtreeHRSL();
  else if(LR==1)readtreeHRSR();

}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::SetInputParam(string ifname){
  param = new ParamMan(ifname.c_str());
  if(param -> SetVal())cout<<"F1TDC parameter setted : really cool acutually"<<endl; 
}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::makehist(){
  string LorR;
  if(LR==0)     LorR="L";
  else if(LR==1)LorR="R";
  for(int i=0;i<16;i++){
    h_s2s0_beta[i]     = new TH1F(Form("h_s2s0_beta%d",i)    , Form("h_s2s0_beta%d",i)     ,200,  -1,1.5);
    h_s2s0_tof[i]      = new TH1F(Form("h_s2s0_tof%d",i)     , Form("h_s2s0_tof%d",i)      ,800,-100,100);
    h_s2s0_tdiff[i]    = new TH1F(Form("h_s2s0_tdiff%d",i)   , Form("h_s2s0_tdiff%d",i)    ,800,-100,100);
    h_s2s0_beta_FB[i]  = new TH1F(Form("h_s2s0_beta_FB%d",i) , Form("h_s2s0_beta_FB%d",i)  ,200,  -1,1.5);
    h_s2s0_tof_FB[i]   = new TH1F(Form("h_s2s0_tof_FB%d",i)  , Form("h_s2s0_tof_FB%d",i)   ,800,-100,100);
    h_s2s0_tdiff_FB[i] = new TH1F(Form("h_s2s0_tdiff_FB%d",i), Form("h_s2s0_tdiff_FB%d",i) ,400, -10, 10);
    set->SetTH1(h_s2s0_beta[i]      ,Form("#beta S2%s%d - S0(F1)"              ,LorR.c_str(),i),"#beta"    ,"counts");
    set->SetTH1(h_s2s0_tof[i]       ,Form("ToF S2%s%d - S0"                    ,LorR.c_str(),i),"ToF[ns]"  ,"counts");
    set->SetTH1(h_s2s0_tdiff[i]     ,Form("TDiff (S2%s%d - S0) - ToF calc"     ,LorR.c_str(),i),"Tdiff[ns]","counts");
    set->SetTH1(h_s2s0_beta_FB[i]   ,Form("#beta S2%s%d - S0 (FB)"             ,LorR.c_str(),i),"#beta"    ,"counts");
    set->SetTH1(h_s2s0_tof_FB[i]    ,Form("ToF S2%s%d - S0 (FB)"               ,LorR.c_str(),i),"ToF[ns]"  ,"counts");
    set->SetTH1(h_s2s0_tdiff_FB[i]  ,Form("TDiff (S2%s%d - S0) - ToF calc (FB)",LorR.c_str(),i),"Tdiff[ns]","counts");

    h2_s2s0_f1tof_fbtof[i]    = new TH2F(Form("h2_s2s0_f1tof_fbtof%d",i)  ,Form("h2_s2s0_f1tof_fbtof%d",i), 100, -20, 20, 100, -20, 20);
    h2_s2s0_f1beta_fbbeta[i]  = new TH2F(Form("h2_s2s0_f1beta_fbbeta%d",i),Form("h2_s2s0_f1beta_fbbeta%d",i), 200, -1, 1.5, 200, -1, 1.5);
    set->SetTH2(h2_s2s0_f1tof_fbtof[i],Form("ToF F1 vs Fbus S2%s%d",LorR.c_str(),i),"ToF(F1)[ns]","ToF(Fbus)[ns]");
    set->SetTH2(h2_s2s0_f1beta_fbbeta[i],Form("#beta F1 vs Fbus S2%s%d",LorR.c_str(),i),"#beta(F1)","#beta(Fbus)");
  }

}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%1000==0)cout<<n <<" / "<<ENum<<endl;
    tree->GetEntry(n);
    if(LR==0 && anaL_oneevent());
    else if(LR==1&&anaR_oneevent());
    else continue;

    for(int i=0;i<16;i++){
      for(int j=0;j<tr_n;j++){
        bool F1Hits = false, FbusHits = false;
        if(s2_trpad[j]==i){
          if(S2T_F1TDC[i]>0 && S2B_F1TDC[i]>0 && S0T_F1TDC[0]>0. && S0B_F1TDC[0]>0.)F1Hits=true;
          if(S2_lt[i]>0. && S2_rt[i]>0. && S0_lt[0]>0 && S0_rt[0]>0)FbusHits=true;
          if(F1Hits){
            h_s2s0_tdiff[i] ->Fill(S2_F1time[i] - S0_F1time[0] - paths2s0[j]/LightVelocity);
            h_s2s0_tof[i]   ->Fill(S2_F1time[i] - S0_F1time[0]);
            h_s2s0_beta[i]  ->Fill(betaF1[j]);
          }
          if(FbusHits){
            h_s2s0_tdiff_FB[i] ->Fill(S2_time[i] - S0_time[0] + paths2s0[j]/LightVelocity);
            h_s2s0_tof_FB[i]   ->Fill(S2_time[i] - S0_time[0]);
            h_s2s0_beta_FB[i]  ->Fill(beta[j]);
          }
          if(F1Hits&&FbusHits){
            h2_s2s0_f1tof_fbtof[i]   ->Fill(S2_F1time[i] - S0_F1time[0],S2_time[i] - S0_time[0]);
            h2_s2s0_f1beta_fbbeta[i] ->Fill(betaF1[j],beta[j]);
          }
        }
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::fit(){

  ////////
  //Fbus//
  ////////
  h_frame[0] = new TH2F("h_frame0","h_frame0",10, -0.8, 15.8,10,-0.5,0.5);
  h_frame[1] = new TH2F("h_frame1","h_frame1",10, -0.8, 15.8,10, 0.0,1.0);
  set->SetTH2(h_frame[0] , "TDiff peak pos each S2(Fbus)","S2 paddle","TDiff peak[ns]" );
  set->SetTH2(h_frame[1] , "TDiff width each S2(Fbus)"   ,"S2 paddle","TDiff width[ns]");
  tg_tdiffFB_pos = new TGraphErrors(); 
  tg_tdiffFB_wid = new TGraphErrors();
  set->SetGrErr(tg_tdiffFB_pos, "TDiff peak pos each S2","S2 paddle","TDiff peak[ns]" ,1,4,23);
  set->SetGrErr(tg_tdiffFB_wid, "TDiff width each S2"   ,"S2 paddle","TDiff width[ns]",1,4,24);

  for(int i=0;i<16;i++){
    ga_tdiffFB[i] = new TF1(Form("ga_tdiffFB%d",i+1),"gaus",-2,2);
    set->SetTF1(ga_tdiffFB[i],2,1,1);
    double min=-50,max=50;
    min = h_s2s0_tdiff_FB[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff_FB[i]->GetMaximumBin()) -3.;
    max = h_s2s0_tdiff_FB[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff_FB[i]->GetMaximumBin()) +3.;
    set->FitGaus(h_s2s0_tdiff_FB[i],min,max,1.5,5);
    h_s2s0_tdiff_FB[i]->Fit(ga_tdiffFB[i],"QR","",min,max);
    tdiffFB_pos[i]  = ga_tdiffFB[i]->GetParameter(1);
    tdiffFB_wid[i]  = ga_tdiffFB[i]->GetParameter(2);
    etdiffFB_pos[i] = ga_tdiffFB[i]->GetParError(1);
    etdiffFB_wid[i] = ga_tdiffFB[i]->GetParError(2);;

    tg_tdiffFB_pos ->SetPoint(i,i,tdiffFB_pos[i]); 
    tg_tdiffFB_wid ->SetPoint(i,i,tdiffFB_wid[i]);
    tg_tdiffFB_pos ->SetPointError(i,0,etdiffFB_pos[i]); 
    tg_tdiffFB_wid ->SetPointError(i,0,etdiffFB_wid[i]);
    
    param->SetTimeTune(CID_FbS2,i,LR,0,tdiffFB_pos[i]);
    param->SetTimeTune(CID_FbS2,i,LR,1,tdiffFB_pos[i]);
  }

  //////
  //F1//
  //////
  tg_tdiffF1_pos = new TGraphErrors(); 
  tg_tdiffF1_wid = new TGraphErrors();
  set->SetGrErr(tg_tdiffF1_pos, "TDiff peak pos each S2(F1)","S2 paddle","TDiff peak[ns]" ,1,4,23);
  set->SetGrErr(tg_tdiffF1_wid, "TDiff width each S2(F1)"   ,"S2 paddle","TDiff width[ns]",1,4,24);
  h_frame[2] = new TH2F("h_frame0","h_frame0",10, -0.8, 15.8,10,-0.5,0.5);
  h_frame[3] = new TH2F("h_frame1","h_frame1",10, -0.8, 15.8,10, 0.0,1.0);
  set->SetTH2(h_frame[2] , "TDiff peak pos each S2(F1TDC)","S2 paddle","TDiff peak[ns]" );
  set->SetTH2(h_frame[3] , "TDiff width each S2(F1TDC)"   ,"S2 paddle","TDiff width[ns]");

  for(int i=0;i<16;i++){
    ga_tdiffF1[i] = new TF1(Form("ga_tdiffF1%d",i+1),"gaus",-2,2);
    set->SetTF1(ga_tdiffF1[i],2,1,1);
    double min=-50,max=50;
    min = h_s2s0_tof[i]->GetXaxis()->GetBinCenter(h_s2s0_tof[i]->GetMaximumBin()) -3.;
    max = h_s2s0_tof[i]->GetXaxis()->GetBinCenter(h_s2s0_tof[i]->GetMaximumBin()) +3.;
    set->FitGaus(h_s2s0_tof[i],min,max,1.5,5);
    h_s2s0_tof[i]->Fit(ga_tdiffF1[i],"QR","",min,max);
    tdiffF1_pos[i]  = ga_tdiffF1[i]->GetParameter(1);
    tdiffF1_wid[i]  = ga_tdiffF1[i]->GetParameter(2);
    etdiffF1_pos[i] = ga_tdiffF1[i]->GetParError(1);
    etdiffF1_wid[i] = ga_tdiffF1[i]->GetParError(2);;

    tg_tdiffF1_pos ->SetPoint(i,i,tdiffF1_pos[i]); 
    tg_tdiffF1_wid ->SetPoint(i,i,tdiffF1_wid[i]);
    tg_tdiffF1_pos ->SetPointError(i,0,etdiffF1_pos[i]); 
    tg_tdiffF1_wid ->SetPointError(i,0,etdiffF1_wid[i]);
  }
  //for(int i=0;i<16;i++){
  //  ga_tdiffF1[i] = new TF1(Form("ga_tdiffF1%d",i+1),"gaus",-2,2);
  //  set->SetTF1(ga_tdiffF1[i],2,1,1);
  //  double min=-50,max=50;
  //  min = h_s2s0_tdiff[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff[i]->GetMaximumBin()) -3.;
  //  max = h_s2s0_tdiff[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff[i]->GetMaximumBin()) +3.;
  //  set->FitGaus(h_s2s0_tdiff[i],min,max,1.5,5);
  //  h_s2s0_tdiff[i]->Fit(ga_tdiffF1[i],"QR","",min,max);
  //  tdiffF1_pos[i]  = ga_tdiffF1[i]->GetParameter(1);
  //  tdiffF1_wid[i]  = ga_tdiffF1[i]->GetParameter(2);
  //  etdiffF1_pos[i] = ga_tdiffF1[i]->GetParError(1);
  //  etdiffF1_wid[i] = ga_tdiffF1[i]->GetParError(2);;

  //  tg_tdiffF1_pos ->SetPoint(i,i,tdiffF1_pos[i]); 
  //  tg_tdiffF1_wid ->SetPoint(i,i,tdiffF1_wid[i]);
  //  tg_tdiffF1_pos ->SetPointError(i,0,etdiffF1_pos[i]); 
  //  tg_tdiffF1_wid ->SetPointError(i,0,etdiffF1_wid[i]);
  //  
  //  param->SetTimeTune(CID_F1S2,i,LR,0,tdiffF1_pos[i]);
  //  param->SetTimeTune(CID_F1S2,i,LR,1,tdiffF1_pos[i]);
  //}
}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::draw(){

  c[0]->Clear();c[0]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[0]->cd(i+1);gPad->SetLogy(1);h_s2s0_tof[i]->Draw();
  }

  c[1]->Clear();c[1]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[1]->cd(i+1);gPad->SetLogy(1);h_s2s0_tdiff[i]->Draw();
  }

  c[2]->Clear();c[2]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[2]->cd(i+1);gPad->SetLogy(1);h_s2s0_tof_FB[i]->Draw();
  }

  c[3]->Clear();c[3]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[3]->cd(i+1);gPad->SetLogy(1);h_s2s0_tdiff_FB[i]->Draw();
  }

  c[4]->Clear();c[4]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[4]->cd(i+1);gPad->SetLogz(1);h2_s2s0_f1tof_fbtof[i] ->Draw("colz");
  }
  c[5]->Clear();c[5]->Divide(2,2);
  c[5]->cd(1);h_frame[0]->Draw("");tg_tdiffFB_pos ->Draw("sameP");///
  c[5]->cd(2);h_frame[1]->Draw("");tg_tdiffFB_wid ->Draw("sameP");///
  c[5]->cd(3);h_frame[2]->Draw("");tg_tdiffF1_pos ->Draw("sameP");///
  c[5]->cd(4);h_frame[3]->Draw("");tg_tdiffF1_wid ->Draw("sameP");///

  c[6]->Clear();c[6]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[6]->cd(i+1);gPad->SetLogy(1);h_s2s0_beta[i]->Draw();
  }

  c[7]->Clear();c[7]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[7]->cd(i+1);gPad->SetLogy(1);h_s2s0_beta_FB[i]->Draw();
  }

  c[8]->Clear();c[8]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[8]->cd(i+1);gPad->SetLogz(1);h2_s2s0_f1beta_fbbeta[i] ->Draw("colz");
  }

  param->WriteToFile("param/tmp.param");
  
}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::savecanvas(string ofname){
  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[NCanvas-1]->Print(Form("%s]",ofname.c_str()) );
  cout<<ofname<<" saved"<<endl;
}
////////////////////////////////////////////////////////////////////////////
bool s2_t0_calib::anaL_oneevent(){
return false;
}
////////////////////////////////////////////////////////////////////////////
bool s2_t0_calib::anaR_oneevent(){
  //cout<<"s2_t0_calib::anaR_oneevent"<<endl;

  convertF1TDCR(param);

  tr_n = (int)R_tr_n;
  //cout<<"tr_n" <<tr_n<<endl;
  if(tr_n>MAX)tr_n=MAX;

  for(int i=0;i<tr_n;i++){
  //cout<<"s2_trpad : "<<R_s2_trpad[i]<<endl;
    s2_trpad[i]=R_s2_trpad[i];
    paths2s0[i]=R_s2_trpath[i] - R_s0_trpath[i];
    beta[i]    =R_tr_beta[i];
    betaF1[i]  =GetBeta_S0S2wF1TDCR(i);
  }

  for(int i=0;i<RS2;i++){
    //F1TDC
    S2T_F1TDC[i] = RS2T_F1TDC[i];
    S2B_F1TDC[i] = RS2B_F1TDC[i];
    S2_F1time[i] = RS2_F1time[i];
    //Fbus TDC
    S2_time[i]   = 1.e+9*R_s2_time[i];
    S2_lt[i]     = R_s2_lt[i];
    S2_rt[i]     = R_s2_rt[i];
  }
  for(int i=0;i<RS0;i++){
    //F1TDC
    S0T_F1TDC[i] = RS0T_F1TDC[i];
    S0B_F1TDC[i] = RS0B_F1TDC[i];
    S0_F1time[i] = RS0_F1time[i];
    //Fbus TDC
    S0_time[i]   = 1.e+9*R_s0_time[i];
    S0_lt[i]     = R_s0_lt[i];
    S0_rt[i]     = R_s0_rt[i];
  }


  if(DR_T4>0. ) return true;
  else return false;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/cosmic1020.root";
  string ofname = "toyamacro/pdf/s2_t0_calib1020.pdf";
  string paramname = "param/default.param";
  int ch;
  int lr=0;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bcop:LR"))!=-1){
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
    case 'L':
      lr = 0;
      break;
    case 'R':
      lr = 1;
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
  s2_t0_calib *calib = new s2_t0_calib();

  calib->SetMaxEvent(MaxNum);
  calib->SetInputParam(paramname);
  calib->SetLR(lr);
  calib->SetRoot(ifname);
  calib->makehist();
  calib->loop();
  calib->fit();
  calib->draw();
  if(output_flag)calib->savecanvas(ofname);
  delete calib;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

