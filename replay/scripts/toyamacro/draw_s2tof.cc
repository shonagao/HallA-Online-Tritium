//////////////////////////////////////////////////
// Check S2 time walk correction parameter(Fbus)//
// by Y. Toyama Oct. 2018                       //
//////////////////////////////////////////////////

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
#include "define.h"

#define Calibration


const int NCanvas = 9;//num of canvas

int main(int argc, char** argv){
  string ifname = "list/run.txt";
  string ofname = "pdf/Rs2_twc.pdf";
  int ch;
  int lr=0;
  int MaxNum = 1000000;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:LR"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input run list: "<<ifname<<endl;
      break;
    case 'w':
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
      cout<<"-f : input runlist filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-L : LHRS"<<endl;
      cout<<"-R : RHRS"<<endl;
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

  string LorR;
  if(lr==0)     LorR="L";
  else if(lr==1)LorR="R";

  TApplication *theApp = new TApplication("App", &argc, argv);
  Tree *tr = new Tree();
  Setting *set=new Setting();
  
  TH2F *h2_tof_adct[16], *h2_tof_adcb[16];
  for(int i=0;i<16;i++){
    h2_tof_adct[i] = new TH2F(Form("h2_tof_adct%d",i),Form("h2_tof_adct%d",i) , 100, -100, 1000, 100, -20, 20);
    h2_tof_adcb[i] = new TH2F(Form("h2_tof_adcb%d",i),Form("h2_tof_adcb%d",i) , 100, -100, 1000, 100, -20, 20);
    set->SetTH2(h2_tof_adct[i]    ,Form("ToF vs S2%s%d top ADC"              ,LorR.c_str(),i) ,"ADC[ch]" ,"ToF[ns]" );
    set->SetTH2(h2_tof_adcb[i]    ,Form("ToF vs S2%s%d bottom ADC"           ,LorR.c_str(),i) ,"ADC[ch]" ,"ToF[ns]" );
  }

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    tr->add_tree(runname);

  }
  tr->pack_tree();

  if(lr==0){//Left HRS
    tr->readtreeS0L();
    tr->readtreeS2L();
  }
  else if(lr==1){//Right HRS
    tr->readtreeS0R();
    tr->readtreeS2R();
  }
  else{cout<<"LR is not defined!"<<endl;return -1;}

  
  int EN=tr->ENum;
  double tof,s0time,s2time,s0at,s0ab,s2at,s2ab;
  double s2tt,s2tb;
  //if( MaxNum>0 && MaxNum<EN) EN = MaxNum;
  for(int n=0;n<EN;n++){
    if(n%5000==0)cout<<n <<" / "<<EN<<"  "<<endl;
    tof=-99.;s0at=-999.;s0ab=-999.;s2at=-999.;s2ab=-999.;
    tr->tree->GetEntry(n);

    for(int s2seg=0;s2seg<16;s2seg++){
      if(lr==0 && (int)tr->L_s0_t_pads[0]==0 && s2seg==(int)tr->L_s2_t_pads[0]){//Left HRL
        s0time = s2ns*tr->R_s0_time[0];
        s2time = s2ns*tr->R_s2_time[s2seg];
        tof    = s2time - s0time;
        s0at   = tr->L_s0_ra_p[0];
        s0ab   = tr->L_s0_la_p[0];
        s2at   = tr->L_s2_ra_p[s2seg];
        s2ab   = tr->L_s2_la_p[s2seg];
        s2tt   = s2ns*tr->L_s2_rt_c[s2seg];
        s2tb   = s2ns*tr->L_s2_lt_c[s2seg];
        h2_tof_adct[s2seg] ->Fill(s2at,tof); 
        h2_tof_adcb[s2seg] ->Fill(s2ab,tof);
      }
      else if(lr==1 && tr->R_s0_rt[0]>0. && tr->R_s0_lt[0]>0. && tr->R_s2_rt[s2seg]>0. && tr->R_s2_lt[s2seg]>0. ){//Right HRL
        //cout<<(int)tr->R_s0_t_pads[i]<<" "<<(int)tr->R_s2_t_pads[i]<<endl;
        s0time = s2ns*tr->R_s0_time[0];
        s2time = s2ns*tr->R_s2_time[s2seg];
        tof    = s2time - s0time;
        s0at   = tr->R_s0_ra_p[0];
        s0ab   = tr->R_s0_la_p[0];
        s2at   = tr->R_s2_ra_p[s2seg];
        s2ab   = tr->R_s2_la_p[s2seg];
        s2tt   = s2ns*tr->R_s2_rt_c[s2seg];
        s2tb   = s2ns*tr->R_s2_lt_c[s2seg];
        h2_tof_adct[s2seg] ->Fill(s2at,tof); 
        h2_tof_adcb[s2seg] ->Fill(s2ab,tof);
      }
    }
  
  }

  TCanvas *c[16];
  
  for(int i=0;i<16;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800,900);
  }
  for(int i=0;i<16;i++){
    c[i]->Clear();
    c[i]->Divide(2,1);
    c[i]->cd(1);gPad->SetLogz(1);h2_tof_adct[i] ->Draw("colz");
    c[i]->cd(2);gPad->SetLogz(1);h2_tof_adcb[i] ->Draw("colz");
  }

  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<16;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[16-1]->Print(Form("%s]",ofname.c_str()) );

  return 0;
}

