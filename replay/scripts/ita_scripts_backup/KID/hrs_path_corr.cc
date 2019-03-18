// Author K. Itabashi Aug. 30th
// HRS nnL experiment missing mass analysis

double corr_R_adc(int i);
double corr_R_x(int i);
double corr_R_th(int i);
double corr_R_alig(int i);
double corr_L_adc(int i);
double corr_L_x(int i);
double corr_L_th(int i);
double corr_L_alig(int i);
double s2f1_off(int i,char* ARM);
const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
const double tdc_time=56.23e-3;//[ns]

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




//=====================================================================//
//============================= Main =================================//
//===================================================================//

int main(int argc, char** argv){
  TApplication *theApp =new TApplication("App",&argc,argv);

  TChain*  T=new TChain("T");
  // T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/tritium_ita111160_111208.root");
  // T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/tritium_ita111200_111210.root");
  // T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/tritium_ita111168_111170.root");
  //  T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_111180.root");

  for(int i=111480;i<111485;i++){
    T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",i));
      for(int j=0;j<5;j++){
       T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d_%d.root",i,j));
    }
  }


 //============= Set Branch Status ==================//

  int max=1000; 
  double RF1[max],LF1[max];
  double Rs0r_ac[max],Rs0l_ac[max],Ls0r_ac[max],Ls0l_ac[max];
  double Rs2r_ac[max],Rs2l_ac[max],Ls2r_ac[max],Ls2l_ac[max];
  double Rs0r_tc[max],Rs0l_tc[max],Ls0r_tc[max],Ls0l_tc[max];
  double Rs2r_tc[max],Rs2l_tc[max],Ls2r_tc[max],Ls2l_tc[max];
  double Ra1t[max],Ra1a[max],Ra1a_p[max],Ra1a_c[max],Ra1sum;
  double Ra2t[max],Ra2a[max],Ra2a_p[max],Ra2a_c[max],Ra2sum;

  double La1t[max],La1a[max],La1a_p[max],La1a_c[max],La1sum;
  double La2t[max],La2a[max],La2a_p[max],La2a_c[max],La2sum;
  double Rp[max],Rpx[max],Rpy[max],Lp[max],Lpx[max],Lpy[max];
  double Rth[max],Rph[max],Rx[max],Rvz[max],Lth[max],Lph[max],Lx[max],Lvz[100];
  double Rbeta[max],Lbeta[max];
  double rs2pathl[max],rs0pathl[max],rtrpathl[max];
  double ls2pathl[max],ls0pathl[max],ltrpathl[max];
  double trigger[100];
  double hallap;
  double Rs2tpads[100],Ls2tpads[100];

  T->SetBranchStatus("*",0);  
  T->SetBranchStatus("HALLA_p",1);
  T->SetBranchAddress("HALLA_p",&hallap); 
  T->SetBranchStatus("DR.evtypebits",1);
  T->SetBranchAddress("DR.evtypebits",trigger); 
 //------ Right Arm -------------//
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
  T->SetBranchStatus("R.s0.ra_c",1); // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.ra_c",Rs0r_ac); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.la_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.la_c",Rs0l_ac); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.ra_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.ra_c",Rs2r_ac);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.la_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.la_c",Rs2l_ac);  // Right arm S2 L-PMT  ADC
  
  T->SetBranchStatus("R.s0.rt_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.rt_c",Rs0r_tc); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.lt_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.lt_c",Rs0l_tc); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.rt_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.rt_c",Rs2r_tc);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.lt_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.lt_c",Rs2l_tc);  // Right arm S2 L-PMT  ADC
  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
  // path length//
  T->SetBranchStatus("R.s2.trpath",1); 
  T->SetBranchAddress("R.s2.trpath",rs2pathl); 
  T->SetBranchStatus("R.s0.trpath",1); 
  T->SetBranchAddress("R.s0.trpath",rs0pathl);
  T->SetBranchStatus("R.tr.pathl",1);  
  T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // (AC1)Aerogel Chrenkov Right ARM ADC //   
 T->SetBranchStatus("R.a1.t",1);
 T->SetBranchStatus("R.a1.a",1);
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchStatus("R.a1.a_p",1);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.t",Ra1t);
 T->SetBranchAddress("R.a1.a",Ra1a);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchAddress("R.a1.a_p",Ra1a_p);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 // (AC2)Aerogel Chrenkov Right ARM ADC //                                   
 T->SetBranchStatus("R.a2.t",1);
 T->SetBranchStatus("R.a2.a",1);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchStatus("R.a2.a_p",1);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.t",Ra2t);
 T->SetBranchAddress("R.a2.a",Ra2a);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 T->SetBranchAddress("R.a2.a_p",Ra2a_p);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c);

 // Target positon information //
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.px",1);
 T->SetBranchAddress("R.tr.px",Rpx);
 T->SetBranchStatus("R.tr.py",1);
 T->SetBranchAddress("R.tr.py",Rpy);
 T->SetBranchStatus("R.tr.ph",1);
 T->SetBranchAddress("R.tr.ph",Rph);
 T->SetBranchStatus("R.tr.th",1);
 T->SetBranchAddress("R.tr.th",Rth);
 T->SetBranchStatus("R.tr.x",1);
 T->SetBranchAddress("R.tr.x",Rx);
 T->SetBranchStatus("R.tr.beta",1);    
 T->SetBranchAddress("R.tr.beta",Rbeta); 
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm ---------------//
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
  T->SetBranchStatus("L.s0.ra_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.ra_c",Ls0r_ac); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.la_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.la_c",Ls0l_ac); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.ra_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.ra_c",Ls2r_ac);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.la_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.la_c",Ls2l_ac);  // Left arm S2 L-PMT  ADC

  T->SetBranchStatus("L.s0.rt_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.rt_c",Ls0r_tc); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.lt_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.lt_c",Ls0l_tc); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.rt_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.rt_c",Ls2r_tc);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.lt_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.lt_c",Ls2l_tc);  // Left arm S2 L-PMT  ADC
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
  // path length//
  T->SetBranchStatus("L.s2.trpath",1); 
  T->SetBranchAddress("L.s2.trpath",ls2pathl); 
  T->SetBranchStatus("L.s0.trpath",1); 
  T->SetBranchAddress("L.s0.trpath",ls0pathl); 
  T->SetBranchStatus("L.tr.pathl",1);   
  T->SetBranchAddress("L.tr.pathl",ltrpathl);

  T->SetBranchStatus("L.tr.beta",1);    
  T->SetBranchAddress("L.tr.beta",Lbeta); 

 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.px",1);
 T->SetBranchAddress("L.tr.px",Lpx);
 T->SetBranchStatus("L.tr.py",1);
 T->SetBranchAddress("L.tr.py",Lpy);
 T->SetBranchStatus("L.tr.ph",1);
 T->SetBranchAddress("L.tr.ph",Lph);
 T->SetBranchStatus("L.tr.th",1);
 T->SetBranchAddress("L.tr.th",Lth);
 T->SetBranchStatus("L.tr.x",1);
 T->SetBranchAddress("L.tr.x",Lx);
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);
 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double min_coin=-260;
 double max_coin=-240;
 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double min_beta=0.95;
 double max_beta=0.97;
 int bin_beta=6000;
 double min_adc=10.0;
 double max_adc=20000.;
 int bin_adc=max_adc-min_adc;


 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",
			 // 100000,-1000,1000);
			 bin_coin,min_coin,max_coin);
 
 double min_rpathl=22.0; double max_rpathl=23.; int bin_rpathl=500;
 // double min_lpathl=22.0; double max_lpathl=23.; int bin_lpathl=500;
 TH2F* hcoin_rpathl=new TH2F("hcoin_rpathl","Coinc time vs R Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
TH2F* hcoin_rpathl_cc=new TH2F("hcoin_rpathl_cc","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_rpathl,max_rpathl,50000,-1000,1000);
 double min_lpathl=22.0; double max_lpathl=23.; int bin_lpathl=500;
 TH2F* hcoin_lpathl=new TH2F("hcoin_lpathl","Coinc time vs L Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
 TH2F* hcoin_lpathl_cc=new TH2F("hcoin_lpathl_cc","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,500,-330,-310);
  
TH2F* hcoin_pathl_c=new TH2F("hcoin_pathl_c","Coinc time vs R Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
  
 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 double mh;
 double m2; 
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rtof[16];
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off; 
 double Rs2_tcorr,Ls2_tcorr;
 int i=8;
 int counts=0;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,n;
 double ac1_adc,ac2_adc;
 double tof_r,tof_l;
 ac1_adc=100.;
 ac2_adc=400.;
 double rpathl,lpathl;
 TLine* lac1=new TLine(0,1.5,ac1_adc,ac1_adc);
 TLine* lac2=new TLine(0,1.5,ac2_adc,ac2_adc);
 double corr_R,corr_L;
 double rpath_corr,lpath_corr,pathl_off;
 
 //----- Cut Parameters ----------//
 double coin_cutmin=-248;
 double coin_cutmax=-244; 
 double rpathl_cutmin=22.2;
 double rpathl_cutmax=22.7;
 double lpathl_cutmin=22.3;
 double lpathl_cutmax=22.8;
 double rbeta_cutmin=0.0;
 double rbeta_cutmax=1.0;
 double lbeta_cutmin=0.9;
 double lbeta_cutmax=1.0;
 double Rvz_cutmin=-0.1;
 double Rvz_cutmax= 0.1;
 double Lvz_cutmin=-0.1;
 double Lvz_cutmax= 0.1;
 double Rx_cutmin= -0.4;
 double Rx_cutmax= 0.4;
 //-------------------------------//
bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig;


for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   // if(trigger[0]==32)cut_trig=true;

 Ls2pads=Ls2tpads[0];
 Rs2pads=Rs2tpads[0]; 
 rpathl=rtrpathl[0]-rs2pathl[0]; // R-HRS path length S2 -RF
 rbeta=pk/Ek; 
 
 lpathl=ltrpathl[0]-ls2pathl[0]; // L-HRS path length S2 -RF
 lbeta=pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 
 Rs2_off=s2f1_off(Rs2pads,"R");
 Ls2_off=s2f1_off(Ls2pads,"L");
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l; //coin time

 //====== Cut condition ========================// 
   cut_trig=true;
   cut_rpathl=false;
   cut_lpathl=false;
   cut_coin=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_trig=false;
   coin_trig=false;
   right_trig=false;

    if(trigger[0]==32)coin_trig=true;
   if(trigger[0]==16)right_trig=true;
   if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
   if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
   if(coin_t<coin_cutmax && coin_t>coin_cutmin)cut_coin=true;
   if(Rvz_cutmin<Rvz[0] && Rvz_cutmax<Rvz[0] && Lvz_cutmin<Lvz[0] && Lvz_cutmax<Lvz[0])cut_vz=true;
   if(Rx_cutmin<Rx[0] && Rx[0]<Rx_cutmax)cut_Rx=true;
 //=======================================//




   //==========================================//
   //========= Fill Hist =====================//
   //========================================//




 if(cut_Rs2 && cut_Ls2 && cut_vz ){
   //hcoin_t->Fill(coin_t);
 hcoin_rpathl->Fill(rpathl,coin_t); 
 hcoin_lpathl->Fill(lpathl,coin_t); 
 }

 }



//====== FitSlicesY =============//
 
  TF1* fcoin_gaus=new TF1("fcoin_gaus[%d]","gaus",min_coin,max_coin);;
  TF1* ffit_r=new TF1("ffit_r","[0]*x+[1]",min_rpathl,max_rpathl);
  TF1* ffit_l=new TF1("ffit_l","[0]*x+[1]",min_lpathl,max_lpathl);

 TObjArray aSlices;  
 hcoin_rpathl->FitSlicesY(fcoin_gaus,0,-1,0,"QRG");
 hcoin_lpathl->FitSlicesY(fcoin_gaus,0,-1,0,"QRG");
 TH1F* hslice_rpathl=(TH1F*)gROOT->FindObject("hcoin_rpathl_1");
 TH1F* hslice_lpathl=(TH1F*)gROOT->FindObject("hcoin_lpathl_1");

 hslice_rpathl->Fit("ffit_r","","");
 hslice_lpathl->Fit("ffit_l","","");

 double rp0,rp1,lp0,lp1;
 rp0=ffit_r->GetParameter(0);
 rp1=ffit_r->GetParameter(1);
 lp0=ffit_l->GetParameter(0);
 lp1=ffit_l->GetParameter(1);

 cout<<"R-HRS Get Parameter p0: "<<rp0<<" p1: "<<rp1<<endl;
 cout<<"L-HRS Get Parameter p0: "<<lp0<<" p1: "<<lp1<<endl; 



 //==== Path Length Correction Hist ==========//

 double min_coin_rpc,max_coin_rpc,min_coin_lpc,max_coin_lpc,bin_coin_rpc,bin_coin_lpc;
 min_coin_rpc=-10.; max_coin_rpc=10.; 
 bin_coin_rpc=(max_coin_rpc-min_coin_rpc)/tdc_time; bin_coin_rpc=(int)bin_coin_rpc;
 min_coin_lpc=-10.; max_coin_lpc=10.; 
 bin_coin_lpc=(max_coin_lpc-min_coin_lpc)/tdc_time; bin_coin_lpc=(int)bin_coin_lpc;

 TH2F* hcoin_rpathl_c=new TH2F("hcoin_rpathl_c","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_lpathl,max_lpathl,bin_coin_rpc,min_coin_rpc,max_coin_rpc); 
TH2F* hcoin_lpathl_c=new TH2F("hcoin_lpathl_c","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,bin_coin_rpc,min_coin_rpc,max_coin_rpc);




for(int k=0;k<evnt;k++){
   T->GetEntry(k);

 Ls2pads=Ls2tpads[0];
 Rs2pads=Rs2tpads[0]; 

 rpathl=rtrpathl[0]-rs2pathl[0]; // R-HRS path length S2 -RF
 rpath_corr=rp0*rpathl+rp1;
 rbeta=pk/Ek; 
 // rpath_corr=rpathl/rbeta/c;

 lpathl=ltrpathl[0]-ls2pathl[0]; // L-HRS path length S2 -RF
 lpath_corr=lp0*lpathl+lp1;
 lbeta=pe_/Ee_; 
 //lpath_corr=lpathl/lbeta/c;


 Rs2_off=s2f1_off(Rs2pads,"R");
 Ls2_off=s2f1_off(Ls2pads,"L");
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l; //coin time
 coin_tc=coin_t-rpath_corr+lpath_corr;

 //====== Cut condition ========================// 
   cut_coin=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;

   if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
   if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
   if(Rvz_cutmin<Rvz[0] && Rvz_cutmax<Rvz[0] && Lvz_cutmin<Lvz[0] && Lvz_cutmax<Lvz[0])cut_vz=true;





   //==========================================//
   //========= Fill Hist =====================//
   //========================================//




 if(cut_Rs2 && cut_Ls2 && cut_vz ){
 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 hcoin_rpathl_c->Fill(rpathl,coin_t-rpath_corr);
 hcoin_lpathl_c->Fill(lpathl,coin_t-lpath_corr);
		      //-lpath_corr);
 hcoin_rpathl_cc->Fill(rpathl,coin_tc);
 hcoin_lpathl_cc->Fill(lpathl,coin_tc);
 
}

 }




  
 TCanvas* ck =new TCanvas("ck","ck");
 ck->Divide(1,2);
 ck->cd(1);
 hcoin_t->Draw();
 ck->cd(2);
 hcoin_tc->Draw();
 

 TCanvas *crpathl=new TCanvas("crpathl","crpathl");
 crpathl->Divide(3,2);
 crpathl->cd(1);
 hcoin_rpathl->Draw();
 crpathl->cd(2);
 hcoin_lpathl->Draw(); 

 crpathl->cd(3);
 hcoin_rpathl->Draw();
 hslice_rpathl->Draw("same");
 crpathl->cd(4);
 hcoin_lpathl->Draw(); 
 hslice_lpathl->Draw("same");
 

 crpathl->cd(5);
 hcoin_rpathl_c->Draw(); 
 crpathl->cd(6);
 hcoin_lpathl_c->Draw();
 /*
 crpathl->cd(5);
 hcoin_lpathl_cc->Draw("colz");
 crpathl->cd(6);
 hcoin_lpathl_cc ->Draw("colz");
 */

 theApp->Run();
 return 0;



}


//==============================================//
//========== Defined Function ==================//
//=============================================//

double corr_L_adc(int i){
  double corr_l_adc;

	  Double_t corr_L_adc_para[14] = {- 1.592e-12, -1.24122e-12, -1.18518e-12, -1.16133e-12, -1.24632e-12, -1.22617e-12, -1.02470e-12, -6.57058e-13, -1.14584e-12, -1.3259e-12, -1.816135e-12, -1.15547e-12,  -1.23475e-12, -1.50406e-12};

	  if(i==0 || i==16)corr_l_adc=0.0;
	  else corr_l_adc=corr_L_adc_para[i-1];

	  return corr_l_adc;

}

double corr_L_x(int i){

  double corr_l_x;
  double cLx=-9.15734e-10;
  Double_t corr_L_x_para[14] = {9.5982e-09, 2.39686e-09, 5.50452e-09, 8.67284e-09, 7.88134e-09, 9.39930e-09,   9.09441e-09, 8.13305e-09, 8.36477e-09, 8.74297e-09, 7.745e-09,  5.94972e-09, 6.22836e-09, 5.52765e-09};

	  if(i==0 || i==16)corr_l_x=0.0;
	  else corr_l_x=corr_L_x_para[i-1]+cLx;

  return corr_l_x;

}


double corr_L_th(int i){

  double corr_l_th;
  double cLth=1.75759e-9;
  Double_t corr_L_th_para[14] = {-5.3783e-08,  - 3.32898e-08, -4.14532e-08, -4.08767e-08, -4.07972e-08, -3.63437e-08,  -3.67840e-08, -3.54952e-08, -3.63706e-08,-3.39145e-08, -3.43925e-08,  -3.05521e-08,-3.07010e-08, -3.79624e-08};

	  if(i==0 || i==16)corr_l_th=0.0;
	  else corr_l_th=corr_L_th_para[i-1]+cLth;

  return corr_l_th;
}

double corr_L_alig(int i){

  double corr_l_alig;
  Double_t corr_L_alig_para[14] = {-5.3783e-08,  - 3.32898e-08, -4.14532e-08, -4.08767e-08, -4.07972e-08, -3.63437e-08,  -3.67840e-08, -3.54952e-08, -3.63706e-08,-3.39145e-08, -3.43925e-08,  -3.05521e-08,-3.07010e-08, -3.79624e-08};

	  if(i==0 || i==16)corr_l_alig=0.0;
	  else corr_l_alig=corr_L_alig_para[i-1];

  return corr_l_alig;
}


double corr_R_x(int i){

  double corr_r_x;
  double corr_R_x_para[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};

 double cRx=4.87486e-11;

	  if(i==0 || i==16)corr_r_x=0.0;
	  else corr_r_x=corr_R_x_para[i-1]+cRx;

  return corr_r_x;
}


double corr_R_th(int i){

  double corr_r_th;
  double corr_R_th_para[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};


 double cRth=-3.06204e-9;

	  if(i==0 || i==16)corr_r_th=0.0;
	  else corr_r_th=corr_R_th_para[i-1]+cRth;

  return corr_r_th;
}

double corr_R_adc(int i){

  double corr_r_adc;
  double corr_R_adc_para[14]={-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};

	  if(i==0 || i==16)corr_r_adc=0.0;
	  else corr_r_adc=corr_R_adc_para[i-1];

	  return corr_r_adc;
}

double corr_R_alig(int i){

  double corr_r_alig;

double corr_R_alig_para[14]={-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};

	  if(i==0 || i==16)corr_r_alig=0.0;
	  else corr_r_alig=corr_R_alig_para[i-1];

  return corr_r_alig;
}

double s2f1_off(int i,char* ARM){

  double Rs2R_off[16]={-8371.15,-8402.64,-8420.95,-8418.19,-8367.51,-8382.20,-8378.34,-8382.16,-8390.86,-8391.47,-8388.92,-8384.14,-8327.55,-8343.96,-8359.88,-8363.72};
  double Rs2L_off[16]={-8483.65,-8477.64,-8495.95,-8493.19,-8517.51,-8494.70,-8528.34,-8494.66,-8503.36,-8466.47,-8501.42,-8496.64,-8515.05,-8493.96,-8472.38,-8476.22};

  double Ls2R_off[16]={-12806.67,-12843.12,-12930.73,-12963.17,-12827.68,-12825.15,-13011.38,-13014.46,-13047.93,-12846.83,-13059.55,-13031.05,-13049.57,-12855.07,-13067.45,-12839.87};
  double Ls2L_off[16]={-12506.67,-12543.12,-12443.23,-12438.17,-12565.18,-12562.65,-12411.38,-12414.46,-12372.93,-12584.33,-12384.55,-12356.05,-12337.07,-12555.07,-12354.95,-12539.87};

  double s2f1_offset;
  if(ARM="R")s2f1_offset= Rs2R_off[i]+Rs2L_off[i];
  if(ARM="L")s2f1_offset=Ls2L_off[i]+Ls2R_off[i];
  return s2f1_offset;

}
