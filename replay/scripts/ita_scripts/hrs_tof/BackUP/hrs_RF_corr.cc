//10/30 Author Itabashi
// RF banch moving macro


// Grobal Function //
 const int chmax=17; // channel of S2 PMT 
 const double tdc_time=56.25e-12;//TDC converse ch->sec [sec/ch]


#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
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
#include "TRandom.h"

class RF_corr
{
public:
  RF_corr(){};
  void arm();
  void ttree();
  void tbranch();
  void makehist();
  void fit();
  void fill();
  void draw();


int nrun,ch;
  TH1F* htof[chmax];
  TH2F* htof_pathl[chmax];
  double pathl[chmax],pathl_s2[chmax],pathl_rf[chmax];
  int evnt;
  bool rarm;
  TChain*  T=new TChain("T");
};
//////////////////////////////////////////////////////////////////////////
void RF_corr::ttree(){
//-------- TTree data input ---------------//

  nrun=94003;
  ch=8;
  T->Add(Form("/w/halla-scifs17exp/triton/itabashi/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",nrun));

cout<<"get tree branch "<<endl;

}
////////////////////////////////////////////////////////////////////////
void RF_corr::arm(){
 if(nrun>90000){rarm=true;}else{rarm=false;}
 cout<<"get HRS arm :"<<rarm<<endl;
}
///////////////////////////////////////////////////////////////////////

void tbranch(){


  int max=10000; 
  double F1[max];
  double s0radc;
  double s0ladc;
  double s2radc[max];
  double s2ladc[max];
  double pathl_s2[max];
  double pathl_rf[chmax];
  int trig;

 //============= Set Branch Status ==================//
  

  T->SetBranchStatus("*",0);  
  if(rarm){
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",F1); 
  T->SetBranchStatus("R.s0.ra_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.ra_c",&s0radc); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.la_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.la_c",&s0ladc); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.ra_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.ra_c",s2radc);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.la_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.la_c",s2ladc);  // Right arm S2 L-PMT  ADC
  T->SetBranchStatus("R.s2.trpath",1);
  T->SetBranchAddress("R.s2.trpath",pathl_s2);
  T->SetBranchStatus("R.tr.pathl",pathl_rf);
  T->SetBranchAddress("R.tr.pathl",pathl_rf);
  }else{
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",F1); 
  T->SetBranchStatus("L.s0.ra_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.ra_c",&s0radc); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.la_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.la_c",&s0ladc); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.ra_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.ra_c",s2radc);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.la_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.la_c",s2ladc);  // Left arm S2 L-PMT  ADC
  T->SetBranchStatus("L.s2.trpath",1);
  T->SetBranchAddress("L.s2.trpath",pathl_s2);
  T->SetBranchStatus("L.tr.pathl",1);
  T->SetBranchAddress("L.tr.pathl",pathl_rf);
  }
cout<<"set tree branch "<<endl;

}
/////////////////////////////////////////////////////////////////////
void RF_corr::fill(){

  evnt=T->GetEntries(); 
  cout<<"GetEntries : "<<evnt<<endl;


}
////////////////////////////////////////////////////////////////////


//=====================================================================//
//============================= Main =================================//
//===================================================================//

int main(int argc, char** argv){
  TApplication *theApp =new TApplication("App",&argc,argv);
  
  RF_corr *HRS_RF_corr=new RF_corr();
  HRS_RF_corr->ttree();
  HRS_RF_corr->arm(); 
  HRS_RF_corr->tbranch();
  HRS_RF_corr->fill();



theApp->Run();
 return 0;
}// END Main
