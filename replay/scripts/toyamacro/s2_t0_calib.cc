///////////////////////////////
//t0 calibration of S2(F1TDC)//
// by Y. Toyama Oct. 2018    //
///////////////////////////////

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

const int NCanvas = 2;//num of canvas

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

    TH1F *h_s2s0_tof[16], *h_s2s0_tdiff[16];
    int LR;//L = 0, R = 1
    int run_num;
    Setting *set;
    ParamMan *param;
    TCanvas *c[NCanvas];


    int tr_n;//num. of track
    double betaF1[MAX],s2_trpad[MAX],paths2s0[MAX];
    double S2_F1time[16],S0_F1time[1];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s2_t0_calib::s2_t0_calib()
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
    h_s2s0_tof[i]   = new TH1F(Form("h_s2s0_tof%d",i)  , Form("h_s2s0_tof%d",i)   ,800,-100,100);
    h_s2s0_tdiff[i] = new TH1F(Form("h_s2s0_tdiff%d",i), Form("h_s2s0_tdiff%d",i) ,800,-100,100);
    set->SetTH1(h_s2s0_tof[i]    ,Form("ToF S2%s%d - S0"               ,LorR.c_str(),i+1),"ToF[ns]","counts");
    set->SetTH1(h_s2s0_tdiff[i]  ,Form("TDiff (S2%s%d - S0) - ToF calc",LorR.c_str(),i+1),"ToF[ns]","counts");
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
        if(s2_trpad[j]==i){
          //cout<<i<<" : "<<S2_F1time[i] <<" ns"<<endl;
          //h_s2s0_tof[i]   ->Fill(S2_F1time[i]);
          h_s2s0_tdiff[i] ->Fill(S0_F1time[0]);
          h_s2s0_tof[i]   ->Fill(S2_F1time[i] - S0_F1time[0]);
          //h_s2s0_tdiff[i] ->Fill(S2_F1time[i] - S0_F1time[0]);
        }
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////
void s2_t0_calib::fit(){

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
    betaF1[i]  =GetBeta_S0S2wF1TDCR(i);
  }

  for(int i=0;i<RS2;i++){
    S2_F1time[i] = RS2_F1time[i];
  }
  for(int i=0;i<RS0;i++){
    S0_F1time[i] = RS0_F1time[i];
  }


  if(DR_T5>0. ) return true;
  else return false;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/cosmic1020.root";
  string ofname = "toyamacro/pdf/s2_t0_calib1020.pdf";
  string paramname = "twlk_param/default.param";
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

  //gSystem->Exit(1);
  theApp->Run();
  return 0;
}

