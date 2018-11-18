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

//============ Opiton =============//
//--- Print ------//
 bool draw=true;
 bool print=false;

 if(draw==0)gROOT->SetBatch(1);


 // TChain* T=new TChain("tree");
 //T->Add("gogami_rootfiles/coin_H2_1.root");  

TChain*  T=new TChain("T");

  for(int i=111160;i<111180;i++){
    T->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/nnL_smallroot/tritium_%d.root",i));
      for(int j=0;j<5;j++){
	T->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/nnL_smallroot/tritium_%d_%d.root",i,j));
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
  double trigger[100],DRT5[10];
  double hallap;
  double Rs2tpads[100],Ls2tpads[100];
  double Rs2trpad[100],Ls2trpad[100];

  T->SetBranchStatus("*",0);  
 //------ Right Arm -------------//

  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
  T->SetBranchStatus("R.s2.trpad",1);
  T->SetBranchAddress("R.s2.trpad",Rs2trpad);
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
  // path length//
  T->SetBranchStatus("R.s2.trpath",1); 
  T->SetBranchAddress("R.s2.trpath",rs2pathl); 
  T->SetBranchStatus("R.tr.pathl",1);  
  T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // Target positon information //
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm ---------------//
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
  T->SetBranchStatus("L.s2.trpad",1);
  T->SetBranchAddress("L.s2.trpad",Ls2trpad);
  // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);
 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double min_coin=-100;
 double max_coin=100.0;
 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 
 double min_coin_c=-100;
 double max_coin_c=100.0;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 

 double min_beta=0.95;
 double max_beta=0.97;
 int bin_beta=6000;

 double min_adc=-500.0;
 double max_adc=20000.;
 int bin_adc=max_adc-min_adc;

 double min_ac1=-500.0;
 double max_ac1=5000.;
 int bin_ac1=max_ac1-min_ac1; 

 double min_ac2=-500.0;
 double max_ac2=20000.;
 int bin_ac2=max_ac2-min_ac2; 


 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t1=new TH1F("hcoin_t1","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t2=new TH1F("hcoin_t2","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t3=new TH1F("hcoin_t3","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);

 double min_rpathl=28.0; double max_rpathl=30.; int bin_rpathl=500;
 double min_lpathl=28.0; double max_lpathl=30.; int bin_lpathl=500; 
  
 double min_coin_rc=80;
 double max_coin_rc=100;
 double bin_coin_rc=(max_coin_rc-min_coin_rc)/tdc_time;
        bin_coin_rc=(int)bin_coin_rc;

 double min_coin_lc=-120;
 double max_coin_lc=-100;
 double bin_coin_lc=(max_coin_lc-min_coin_lc)/tdc_time;
        bin_coin_lc=(int)bin_coin_lc;
 
 TH2F* hcoin_rpathl=new TH2F("hcoin_rpathl","Coinc time vs R Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
 TH2F* hcoin_lpathl=new TH2F("hcoin_lpathl","Coinc time vs L Path Length Hist ",bin_lpathl,min_lpathl,max_lpathl,bin_coin,min_coin,max_coin); 

 TH2F* hcoin_rpathl_c=new TH2F("hcoin_rpathl_c","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_rpathl,max_rpathl,bin_coin_rc,min_coin_rc,max_coin_rc); 
 TH2F* hcoin_lpathl_c=new TH2F("hcoin_lpathl_c","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,bin_coin_lc,min_coin_lc,max_coin_lc);

 TH2F* hcoin_rpathl_cc=new TH2F("hcoin_rpathl_cc","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_rpathl,max_rpathl,bin_coin_c,min_coin_c,max_coin_c);
 TH2F* hcoin_lpathl_cc=new TH2F("hcoin_lpathl_cc","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,bin_coin_c,min_coin_c,max_coin_c);
  
 TH2F* hcoin_ac1=new TH2F("hcoin_ac1","Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 TH2F* hcoin_ac2=new TH2F("hcoin_ac2","Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);

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
 int nac1,nac2,nac3,n;
 double ac1_adc,ac2_adc,ac2_adc2;
 double tof_r,tof_l;
 ac1_adc=150.;
 ac2_adc=5000.;
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr,pathl_off;
 double mm;
 //----- Cut Parameters ----------//
 double coin_cutmin=-248;
 double coin_cutmax=-244; 
 double rpathl_cutmin=28.7;
 double rpathl_cutmax=29.4;
 double lpathl_cutmin=28.6;
 double lpathl_cutmax=29.2;
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
 bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig,cut_track,cut_s0;


for(int k=0;k<evnt;k++){
   T->GetEntry(k);

 

 pe_=Lp[0];//*sqrt(1+pow(Lth[0],2)+pow(Lph[0],2));
 pk=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 ppi=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 pe=4.313; // GeV   //hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Epi=sqrt(pow(ppi,2)+pow(mpi,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 Ls2pads=Ls2tpads[0];
 Rs2pads=Rs2tpads[0];
 

 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rbeta=pk/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 
 Rs2_off=s2f1_off(Rs2pads,"R");
 Ls2_off=s2f1_off(Ls2pads,"L");
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l; //coin time
 coin_tc=coin_t+rpath_corr-lpath_corr;//+pathl_off;

//====== Cut condition ========================// 
   cut_ac1=false;
   cut_ac2=false;
   // cut_trig=true;
   // cut_beta=false;
   cut_rpathl=false;
   cut_lpathl=false;
   // cut_coin=false;
   // cut_rbeta=false;
   // cut_lbeta=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   // cut_Rx=false;
   // cut_trig=false;
   // coin_trig=true;
   // right_trig=false;
   cut_track=false;
   cut_s0=false;
   if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
   if(Ra1sum<ac1_adc)cut_ac1=true;
   if(Ra2sum<ac2_adc)cut_ac2=true;
   //   if(rbeta>0.963 && rbeta<0.966)cut_beta=true;
   if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
   if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
   // if(coin_t<coin_cutmax && coin_t>coin_cutmin)cut_coin=true;
   if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
   if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
   // if(rbeta_cutmin<rbeta && rbeta_cutmax)cut_rbeta=true;
   // if(lbeta_cutmin<rbeta && lbeta_cutmax)cut_lbeta=true;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
   // if(Rx_cutmin<Rx[0] && Rx[0]<Rx_cutmax)cut_Rx=true;
   if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
 //=======================================//



   //==========================================//
   //========= Fill Hist =====================//
   //========================================//



 if(cut_Rs2 && cut_Ls2 && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){

 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 // hmm->Fill(mh);
 hcoin_rpathl->Fill(rpathl,coin_t);
 hcoin_rpathl_c->Fill(rpathl,coin_t+rpath_corr);
 hcoin_rpathl_cc->Fill(rpathl,coin_tc);
 hcoin_lpathl->Fill(lpathl,coin_t);
 hcoin_lpathl_c->Fill(lpathl,coin_t-lpath_corr);
 hcoin_lpathl_cc->Fill(lpathl,coin_tc);
 hcoin_ac1->Fill(coin_tc,Ra1sum);
 hcoin_ac2->Fill(coin_tc,Ra2sum);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2
 }
 if(cut_Rs2 && cut_Ls2 && cut_ac1 && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t1->Fill(coin_tc); // AC1 cut

 }
 if(cut_Rs2 && cut_Ls2 && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t2->Fill(coin_tc); //AC2 cut

 }
 if(cut_Rs2 && cut_Ls2  && cut_ac1 && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t3->Fill(coin_tc); //AC1 & AC2 cut
 
  }
 
 }




 n=hcoin_t->GetEntries();
 nac1=hcoin_t1->GetEntries();
 nac2=hcoin_t2->GetEntries();
 nac3=hcoin_t3->GetEntries();
  
 TLine* lac=new TLine(min_coin_c,ac1_adc,max_coin_c,ac1_adc);
 lac->SetLineWidth(2);
 lac->SetLineColor(2);
 
 TCanvas* ccoin_ac=new TCanvas("ccoin_ac","ccoin_ac");
 ccoin_ac->Divide(1,2);
 ccoin_ac->cd(1);
 hcoin_ac1->Draw("colz");
 lac->DrawLine(min_coin_c,ac1_adc,max_coin_c,ac1_adc);
 ccoin_ac->cd(2);
 hcoin_ac2->Draw("colz");
 lac->DrawLine(min_coin_c,ac2_adc,max_coin_c,ac2_adc); 
 TCanvas* ccoin=new TCanvas("ccoin","ccoin");
 ccoin->cd();
 hcoin_tc->SetLineColor(6);
 hcoin_tc->Draw();
 hcoin_t1->SetLineColor(2);
 hcoin_t2->SetLineColor(3);
 hcoin_t3->SetLineColor(1);
 hcoin_t1->Draw("same");
 hcoin_t2->Draw("same");
 hcoin_t3->Draw("same");

 //============== Efficiency analysis =======================//

 /*
 int bin_ac1_adc=hcoin_ac1->GetXaxis()->FindBin(ac1_adc);
 int bin_max_ac1=hcoin_ac1->GetXaxis()->FindBin(max_ac1);
 TH1F* hcoin_ac1_p =hcoin_ac1->ProjectionX("hcoin_ac1_p",bin_ac1_adc,bin_max_adc);
 int bin_ac2_adc=hcoin_ac2->GetXaxis()->FindBin(ac2_adc);
 int bin_max_ac2=hcoin_ac2->GetXaxis()->FindBin(max_ac2);
 TH1F* hcoin_ac2_p =hcoin_ac2->ProjectionX("hcoin_ac2_p",bin_ac2_adc,bin_max_adc);
 */

 TH1F* hcoin_tx=(TH1F*)hcoin_t2->Clone();
 hcoin_tx->SetName("hcoin_tx");


 TF1* facc=new TF1("facc","[0]",min_coin_c,max_coin_c);
 hcoin_tx->Fit("facc","R","",min_coin_c,-18.);
 double p0_acc=facc->GetParameter(0);
 // double p1_acc=facc->GetParameter(1);
 TF1* fp=new TF1("fp","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp->FixParameter(3,p0_acc);
 fp->SetParameter(1,-14);
 fp->SetParameter(2,0.788);
 hcoin_tx->Fit("fp","R","",-16.,-12.);
 double n_p=fp->GetParameter(0);
 double mean_p=fp->GetParameter(1);
 double sig_p=fp->GetParameter(2);
 TF1* fpi=new TF1("fpi","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi->SetParameter(3,p0_acc);
 fpi->SetParameter(1,-2.59);
 fpi->SetParameter(2,0.384);
 hcoin_tx->Fit("fpi","R","",-4.,-1.);
 double n_pi=fpi->GetParameter(0);
 double mean_pi=fpi->GetParameter(1);
 double sig_pi=fpi->GetParameter(2);
 TF1* fk=new TF1("fk","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk->SetParameter(3,p0_acc);
 fk->SetParameter(1,-5.5);
 fk->SetParameter(2,0.644);
 fk->SetParLimits(2,0,0.65);
 hcoin_tx->Fit("fk","R","",-7.,-4.5);
 double n_k=fk->GetParameter(0);
 double mean_k=fk->GetParameter(1);
 double sig_k=fk->GetParameter(2);


 TF1* fcoin=new TF1("fcoin","gausn(0)+gausn(3)+gausn(6)+pol0(9)");
 fcoin->SetParameters(n_pi,mean_pi,sig_pi,n_k,mean_k,sig_k,n_p,mean_p,sig_p,p0_acc);
 fcoin->SetParLimits(4,mean_k-0.5*sig_k,mean_k+0.5*sig_k);
 fcoin->SetParLimits(5,0.8*sig_k,1.2*sig_k);
 hcoin_tx->Fit("fcoin","","",min_coin_c,max_coin_c);

n_pi=fcoin->GetParameter(0); mean_pi=fcoin->GetParameter(1);  sig_pi=fcoin->GetParameter(2);
n_k=fcoin->GetParameter(3); mean_k=fcoin->GetParameter(4);  sig_k=fcoin->GetParameter(5);
n_p=fcoin->GetParameter(6); mean_p=fcoin->GetParameter(7);  sig_p=fcoin->GetParameter(8);
 p0_acc=fcoin->GetParameter(9);
 // p1_acc=fcoin->GetParameter(10);


 double sum_k=n_k/tdc_time; 
 double sum_pi=n_pi/tdc_time;
 double sum_p=n_p/tdc_time;

 fk->SetParameters(n_k,mean_k,sig_k);
 fp->SetParameters(n_p,mean_p,sig_p);
 fpi->SetParameters(n_pi,mean_pi,sig_pi);

 int bin_pi_min=hcoin_ac2->GetXaxis()->FindBin(mean_pi-3*sig_pi);
 int bin_pi_max=hcoin_ac2->GetXaxis()->FindBin(mean_pi+3*sig_pi);
 int bin_k_min=hcoin_ac2->GetXaxis()->FindBin(-5.55-3*0.644);
 int bin_k_max=hcoin_ac2->GetXaxis()->FindBin(-5.55+3*0.644);
 int bin_coint_min=hcoin_ac2->GetXaxis()->FindBin(min_coin_c);
 int bin_coint_max=hcoin_ac2->GetXaxis()->FindBin(max_coin_c);
 TH1D* hac2_pi=hcoin_ac2->ProjectionY("hac2_pi",bin_pi_min,bin_pi_max);
 TH1D* hac2_k=hcoin_ac2->ProjectionY("hac2_k",bin_k_min,bin_k_max);
 TH1D* hac2_all=hcoin_ac2->ProjectionY("hac2_all",bin_coint_min,bin_coint_max);


 //======== Draw TCanvas ==============//

 TCanvas* cAC2_test=new TCanvas("cAC2_test","cAC2_test");
 cAC2_test->cd();
 hac2_all->SetLineColor(1);
 hac2_all->Draw();
 hac2_pi->SetLineColor(46);
 // hac2_pi->SetFillColor(46);
 hac2_pi->Draw("same");
 hac2_k->SetLineColor(4);
 // hac2_k->SetFillColor(4);
 hac2_k->Draw("same");
 TCanvas* c0=new TCanvas("c0","c0");
 c0->cd();
 hcoin_tx->Draw();
 fcoin->SetLineColor(kRed);
 fk->SetLineColor(4);
 fk->SetFillColor(4);
 fk->SetFillStyle(3001);
 fpi->SetFillStyle(3001);
 fpi->SetLineColor(46);
 fpi->SetFillColor(46);
 fp->SetFillStyle(3001);
 fp->SetLineColor(8);
 fp->SetFillColor(8);
 fcoin->Draw("same");
 fk->Draw("same");
 fp->Draw("same");
 fpi->Draw("same");
 TCanvas* cAC=new TCanvas("cAC","cAC");
 cAC->cd();
 ha1_a2->Draw("colz");
 lac->DrawLine(min_ac1,ac2_adc,max_ac1,ac2_adc);
 lac->DrawLine(ac1_adc,min_ac2,ac1_adc,max_ac2);

 TCanvas* c_pathc =new TCanvas("c_pathc","c_pathc");
 c_pathc->Divide(1,2);
 c_pathc->cd(1);
 hcoin_t->Draw();
 c_pathc->cd(2);
 hcoin_tc->Draw();

 TCanvas *crpathl=new TCanvas("crpathl","crpathl");
 crpathl->Divide(2,3);
 crpathl->cd(1);
 hcoin_rpathl->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);
 crpathl->cd(2);
 hcoin_lpathl->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 
 crpathl->cd(3);
 hcoin_rpathl_c->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);
 crpathl->cd(4);
 hcoin_lpathl_c->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 
 crpathl->cd(5);
 hcoin_rpathl_cc->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);
 crpathl->cd(6);
 hcoin_lpathl_cc ->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 

 //================ Print Canvas =================================//
TString name;
 if(print){
   // name.Form("./pdf/hdrogen_run.pdf");
 name.Form("./pdf/hdrogen_test.pdf");
 ccoin->Print(name+"[","pdf");
 ccoin->Print(name,"pdf");
 c0->Print(name,"pdf");
 ccoin_ac->Print(name,"pdf");
 cAC->Print(name,"pdf");
 c_pathc->Print(name,"pdf");
 crpathl->Print(name,"pdf");
 crpathl->Print(name +"]","pdf");
}


 //====== Comment Out =====================//

 cout<<" Event w/o ac cut : "<<n<<endl;
 cout<<" Event w/ ac1 cut : "<<nac1<<endl;
 cout<<" Event w/ ac2 cut : "<<nac2<<endl;
 cout<<" Event w/ ac1 & ac2 cut : "<<nac3<<endl;
 cout<<"======= Get Fit Parameters ============"<<endl;
 cout<<"Accidental BG p0: "<<p0_acc<<endl;
 cout<<"Proton Fit parameters "<<endl;
 cout<<"mean "<<mean_p<<endl;
 cout<<"sigma "<<sig_p<<endl;
 cout<<"Pion Fit Parameters "<<endl;
 cout<<"mean "<<mean_pi<<endl;
 cout<<"sigma "<<sig_pi<<endl;
 cout<<"Kaon Fit Paramters "<<endl;
 cout<<"mean "<<mean_k<<endl;
 cout<<"sigma "<<sig_k<<endl;
 cout<<"Sum of Kaon "<<sum_k<<endl;
 cout<<"Sum of Pion "<<sum_pi<<endl;
 cout<<"Sum of Proton "<<sum_p<<endl;
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
  if(ARM=="R")s2f1_offset=Rs2R_off[i]+Rs2L_off[i];
  else  if(ARM=="L")s2f1_offset=Ls2L_off[i]+Ls2R_off[i];
  else {cout<<"false read out !!"<<endl;}
  return s2f1_offset;

}
