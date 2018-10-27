//////////////////////////////////////
//time walk calibration of S2(Fbus) //
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

#include "Setting.h"

#define Calibration

const int NCanvas = 3;//num of canvas

class s2_twc_calib
{
 public:
         s2_twc_calib();
        ~s2_twc_calib();
  void makehist();
  void loop();
  void fit();
  void search_best();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void SetInputParam(string ifname);
  void SetLR(int lr){LR=lr;}
  int ENum;
  int s2seg;
  double tof,s0time,s2time,s0at,s0ab,s2at,s2ab;

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;

    TH2F *h2_tof_s0at, *h2_tof_s0ab, *h2_tof_s2at, *h2_tof_s2ab;
    TH2F *h2_factor;
    TH1D *h_s2at, *h_s2ab;
    TH2F *h_frame[4];
    int LR;//L = 0, R = 1
    Setting *set;
    TCanvas *c[NCanvas];
    TGraphErrors *tg_tdiffFB_pos, *tg_tdiffFB_wid;
    TGraphErrors *tg_tdiffF1_pos, *tg_tdiffF1_wid;

    TF1 *f_T,*f_B;

    TFile *ifp;
    TTree *tree;

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s2_twc_calib::s2_twc_calib()
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
s2_twc_calib::~s2_twc_calib(){
}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::SetRoot(string ifname){
  ifp = new TFile(ifname.c_str());
  tree = (TTree*)ifp->Get("tree");
  tree->SetBranchAddress("s2seg"  ,&s2seg );
  tree->SetBranchAddress("tof"    ,&tof   );
  tree->SetBranchAddress("s0time" ,&s0time);
  tree->SetBranchAddress("s2time" ,&s2time);
  tree->SetBranchAddress("s0at"   ,&s0at  );
  tree->SetBranchAddress("s0ab"   ,&s0ab  );
  tree->SetBranchAddress("s2at"   ,&s2at  );
  tree->SetBranchAddress("s2ab"   ,&s2ab  );
  ENum = tree->GetEntries();
}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::makehist(){
  string LorR;
  if(LR==0)     LorR="L";
  else if(LR==1)LorR="R";
  h2_tof_s0at  = new TH2F("h2_tof_s0at"  ,"h2_tof_s0at", 100, -100, 4000, 100, -20, 20);
  h2_tof_s0ab  = new TH2F("h2_tof_s0ab"  ,"h2_tof_s0ab", 100, -100, 8000, 100, -20, 20);
  h2_tof_s2at  = new TH2F("h2_tof_s2at"  ,"h2_tof_s2at", 100, -100, 1000, 100, -20, 20);
  h2_tof_s2ab  = new TH2F("h2_tof_s2ab"  ,"h2_tof_s2ab", 100, -100, 1000, 100, -20, 20);
  h2_factor    = new TH2F("h2_factor"    ,"h2_factor"  , 100, 0.4,1.8, 100, 0.4,1.8);

  set->SetTH2(h2_tof_s0at  ,Form("ToF vs S0%s top ADC"        ,LorR.c_str())       ,"ADC[ch]" ,"ToF[ns]" );
  set->SetTH2(h2_tof_s0ab  ,Form("ToF vs S0%s bottom ADC"     ,LorR.c_str())       ,"ADC[ch]" ,"ToF[ns]" );
  set->SetTH2(h2_tof_s2at  ,Form("ToF vs S2%s%d top ADC"      ,LorR.c_str(),s2seg) ,"ADC[ch]" ,"ToF[ns]" );
  set->SetTH2(h2_tof_s2ab  ,Form("ToF vs S2%s%d bottom ADC"   ,LorR.c_str(),s2seg) ,"ADC[ch]" ,"ToF[ns]" );
  set->SetTH2(h2_factor    ,Form("S2%s%d time walk correction",LorR.c_str(),s2seg) ,"Factor T","Factor B");

  f_T = new TF1("f_T","[0]/sqrt(x)-[1]",-100,1000);
  f_B = new TF1("f_B","[0]/sqrt(x)-[1]",-100,1000);
  set->SetTF1(f_T,2,1,1);
  set->SetTF1(f_B,2,1,1);

}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%5000==0)cout<<n <<" / "<<ENum<<endl;
    tree->GetEntry(n);

    if(s0time>205){//to cut double peak(run94010) should be checked and removed if you can
      h2_tof_s0at  ->Fill(s0at,tof);
      h2_tof_s0ab  ->Fill(s0ab,tof);
      h2_tof_s2at  ->Fill(s2at,tof);
      h2_tof_s2ab  ->Fill(s2ab,tof);
    }

  }

}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::fit(){
  f_T->SetParameter(1,0);
  f_B->SetParameter(1,0);
  TF1 *f_ga = new TF1("f_ga","gaus",-20,20);    
  h2_tof_s2at  ->FitSlicesY(f_ga,0,-1,0,"QRG2");
  h_s2at = (TH1D*)gROOT->FindObject("h2_tof_s2at_1");
  h_s2at ->Fit(f_T,"QR","",0,800);

  h2_tof_s2ab  ->FitSlicesY(f_ga,0,-1,0,"QNRG2");
  h_s2ab = (TH1D*)gROOT->FindObject("h2_tof_s2ab_1");
  h_s2ab ->Fit(f_B,"QR","",0,800);

}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::search_best(){
  cout<<"search best parameter"<<endl;

}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::draw(){

  c[0]->Clear();c[0]->Divide(2,2);
  c[0]->cd(1);gPad->SetLogz(1);h2_tof_s0at->Draw("colz");
  c[0]->cd(2);gPad->SetLogz(1);h2_tof_s0ab->Draw("colz");
  c[0]->cd(3);gPad->SetLogz(1);h2_tof_s2at->Draw("colz");
  c[0]->cd(4);gPad->SetLogz(1);h2_tof_s2ab->Draw("colz");

  c[1]->Clear();c[1]->Divide(2,2);
  c[1]->cd(1);gPad->SetLogz(1);h2_tof_s2at->Draw("colz");h_s2at->Draw("samePE");
  c[1]->cd(2);gPad->SetLogz(1);h2_tof_s2ab->Draw("colz");h_s2ab->Draw("samePE");
  c[1]->cd(3);gPad->SetLogz(1);h_s2at->Draw("PE");
  c[1]->cd(4);gPad->SetLogz(1);h_s2ab->Draw("PE");

  c[2]->Clear();c[2]->Divide(1,1);
  c[2]->cd(1);gPad->SetLogz(0);h2_factor->Draw("colz");

}
////////////////////////////////////////////////////////////////////////////
void s2_twc_calib::savecanvas(string ofname){
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
  string ofname = "toyamacro/pdf/s2_twc_calib1020.pdf";
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
  s2_twc_calib *calib = new s2_twc_calib();

  calib->SetMaxEvent(MaxNum);
  calib->SetLR(lr);
  calib->SetRoot(ifname);
  calib->makehist();
  calib->loop();
  calib->fit();
  calib->search_best();
  calib->draw();
  if(output_flag)calib->savecanvas(ofname);
  delete calib;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

