// Author K. Itabashi Nov. 03     //
// R-HRS AC Efficiency study code //

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
  //  T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/tritium_ita111160_111208.root");
  // T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/ita_Rootfiles/tritium_ita111168_111170.root");
  for(int i=111160;i<111210;i++){
    T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",i));
    for(int j=0;j<5;j++){
      T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d_%d.root",i,j));
}
  }

 //============= Set Branch Status ==================//

  int max=100; 
  double RF1[max],LF1[max];
  double Rs0r_ac[max],Rs0l_ac[max],Ls0r_ac[max],Ls0l_ac[max];
  double Rs2r_ac[max],Rs2l_ac[max],Ls2r_ac[max],Ls2l_ac[max];
  double Rs0r_tc[max],Rs0l_tc[max],Ls0r_tc[max],Ls0l_tc[max];
  double Rs2r_tc[max],Rs2l_tc[max],Ls2r_tc[max],Ls2l_tc[max];
  double Ra1t[max],Ra1a[max],Ra1a_p[max],Ra1a_c[max],Ra1sum;
  double Ra2t[max],Ra2a[max],Ra2a_p[max],Ra2a_c[max],Ra2sum;

  double La1t[max],La1a[max],La1a_p[max],La1a_c[max],La1sum;
  double La2t[max],La2a[max],La2a_p[max],La2a_c[max],La2sum;
  double Rp[max],Rpx[max],Rpy[max],Rvz[max],Lp[max],Lpx[max],Lpy[max],Lvz[max];
  double Rth[max],Rph[max],Rx[max],Lth[max],Lph[max],Lx[max];
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


 //========= Hist Definition ============//
 double min_coin=-265;
 double max_coin=-240;
 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double min_coin_pc=-500;
 double max_coin_pc=500;
 int bin_coin_pc=800;
 double min_beta=0.962;
 double max_beta=0.97;
 int bin_beta=6000;
 double min_adc=10.0;
 double max_adc=20000.;
 int bin_adc=max_adc-min_adc;
 double min_ac1=0.0; double max_ac1=5000.;
 int bin_ac1=(max_ac1-min_ac1);
 double min_ac2=0.0; double max_ac2=20000.;
 int bin_ac2=(int)(max_ac2-min_ac2);
 double min_rpathl=22.0; double max_rpathl=23.0;
 double bin_rpathl=(max_rpathl-min_rpathl)/tdc_time/c;
 bin_rpathl=(int)bin_rpathl;
 double min_lpathl=22.0; double max_lpathl=23.0;
 double bin_lpathl=(max_lpathl-min_lpathl)/tdc_time/c;
 bin_lpathl=(int)bin_lpathl;

 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_cut=new TH1F("hcoin_cut","Coincidence time with Pion Select cut S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH2F* hbeta_ac1=new TH2F("hbeta_ac1","R-HRS beta vs AC1 ADC Hist",bin_beta,min_beta,max_beta,bin_ac1,min_ac1,max_ac1); 
 TH2F* hbeta_ac2=new TH2F("hbeta_ac2","R-HRS beta vs AC2 ADC Hist",bin_beta,min_beta,max_beta,bin_ac2,min_ac2,max_ac2); 
 TH2F* hac1_ac2=new TH2F("hac1_ac2","R-HRS AC1 vs AC2 ADC Hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2); 
 TH2F* hcoin_ac1=new TH2F("hcoin_ac1","Coincidence vs AC1 ADC sum Hist",bin_ac1,min_ac1,max_ac1,bin_coin,min_coin,max_coin);
 TH2F* hcoin_ac2=new TH2F("hcoin_ac2","Coincidence vs AC2 ADC sum Hist",bin_ac2,min_ac2,max_ac2,bin_coin,min_coin,max_coin);
 TH2F* hbeta_ac1_cut=new TH2F("hbeta_ac1_cut","R-HRS beta vs AC1 ADC Hist w/ Select pion",bin_beta,min_beta,max_beta,bin_ac1,min_ac1,max_ac1); 
  TH2F* hcoin_rbeta=new TH2F("hcoin_rbeta","Coin-Time vs R beta Hist",bin_beta,min_beta,max_beta,bin_coin,min_coin,max_coin);
 //TH2F* hcoin_ac1=new TH2F("hcoin_ac1","R-HRS Coin-Time vs AC1 ADC Hist",bin_coin,min_coin,max_coin,bin_ac1,min_ac1,max_ac1); 
 //TH2F* hcoin_ac2=new TH2F("hcoin_ac2","R-HRS Coin-Time vs AC2 ADC Hist",bin_coin,min_coin,max_coin,bin_ac2,min_ac2,max_ac2); 
 


 TH2F* hbeta_ac2_cut=new TH2F("hbeta_ac2_cut","R-HRS beta vs AC2 ADC Hist w/ Select pion",bin_beta,min_beta,max_beta,bin_ac2,min_ac2,max_ac2); 
 TH2F* hac1_ac2_cut=new TH2F("hac1_ac2_cut","R-HRS AC1 vs AC2 ADC Hist w/ Select pion",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2); 
 TH2F* hcoin_rpathl=new TH2F("hcoin_rpathl","RHRS Path Length vs Coincidence time Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin);
 TH2F* hcoin_lpathl=new TH2F("hcoin_lpathl","LHRS Path Length vs Coincidence time Hist ",bin_lpathl,min_lpathl,max_lpathl,bin_coin,min_coin,max_coin);
TH2F* hcoin_rpathlc=new TH2F("hcoin_rpathlc","RHRS Path Length vs Coincidence time Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin_pc,min_coin_pc,max_coin_pc);
 TH2F* hcoin_lpathlc=new TH2F("hcoin_lpathlc","LHRS Path Length vs Coincidence time Hist ",bin_lpathl,min_lpathl,max_lpathl,bin_coin_pc,min_coin_pc,max_coin_pc);
 TH1F* hcoin_pc=new TH1F("hcoin_pc","Coincidence tim w/ Path Length Corrction",bin_coin_pc,min_coin_pc,max_coin_pc);


 //======= Define Functions =========//
 double pe_,pe,pk,Ee,Ee_,Epi,Ek;
 int Ls2pads,Rs2pads;
 double coin_t,coin_tc,coin_cr,coin_cl,tof_r,tof_l,lbeta,rbeta,rpathl,lpathl,pathl_off;
 double rpath_corr,lpath_corr;
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
 pathl_off=94;
 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;


 for(int k=0;k<evnt;k++){
 T->GetEntry(k);

 //------ Get Parameters ------------------//
 Rs2pads=(int)Rs2tpads[0];
 Ls2pads=(int)Ls2tpads[0];
 pe_=Lp[0]*sqrt(1-pow(Lth[0],2)+pow(Lph[0],2));
 pk=Rp[0]*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 pe=hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 rpath_corr=rpathl/rbeta/c;
 lpath_corr=lpathl/lbeta/c;
 rpathl=rtrpathl[0]-rs2pathl[0]; // R-HRS path length S2 -RF
 rbeta=pk/Ek; 
 lpathl=ltrpathl[0]-ls2pathl[0]; // L-HRS path length S2 -RF
 lbeta=pe/Ee; 
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37])/2.0)*tdc_time;
 coin_t=tof_r-tof_l; //coin time
 coin_tc=coin_t+rpath_corr+lpath_corr+pathl_off;
 coin_cr=coin_t+rpath_corr;
 coin_cl=coin_t+lpath_corr;
 //-----------------------------------------------//

 //====== Cut condition ========================//
 cut_rpathl=false;
 cut_lpathl=false;
 cut_coin=false;
 cut_rbeta=false;
 cut_lbeta=false;
 cut_Rs2=false;
 cut_Ls2=false;
 cut_vz=false;
 cut_Rx=false;
 cut_trig=false;
 coin_trig=false;
 right_trig=false;
 if(trigger[0]==32)coin_trig=true;
 if(trigger[0]==16)right_trig=true;
 if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
 if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
 if(coin_t<coin_cutmax && coin_t>coin_cutmin)cut_coin=true;
 if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
 if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
 if(rbeta_cutmin<rbeta && rbeta_cutmax)cut_rbeta=true;
 if(lbeta_cutmin<rbeta && lbeta_cutmax)cut_lbeta=true;
 if(Rvz_cutmin<Rvz[0] && Rvz_cutmax<Rvz[0] && Lvz_cutmin<Lvz[0] && Lvz_cutmax<Lvz[0])cut_vz=true;
 if(Rx_cutmin<Rx[0] && Rx[0]<Rx_cutmax)cut_Rx=true;
 //=======================================//


 //============= Fill Hist ===============//

 if(cut_Rs2 && cut_Ls2 && coin_trig){
   //--- Coin Time Hist ----//
 hcoin_t->Fill(coin_t);
 hcoin_rpathl->Fill(rpathl,coin_t);
 hcoin_lpathl->Fill(lpathl,coin_t);
 hcoin_rpathlc->Fill(rpathl,coin_cr);
 hcoin_lpathlc->Fill(lpathl,coin_cl);
 if(cut_rpathl && cut_lpathl && cut_rbeta && cut_lbeta && cut_vz && cut_Rx)hcoin_cut->Fill(coin_t);
 
}
 if(cut_Rs2 && right_trig){
 //----- Beta vs AC ---------//
 hbeta_ac1->Fill(rbeta,Ra1sum);
 hbeta_ac2->Fill(rbeta,Ra2sum);
 hac1_ac2->Fill(Ra1sum,Ra2sum);
 //---- Coin vs AC --------//
 hcoin_ac1->Fill(Ra1sum,coin_t);
 hcoin_ac2->Fill(Ra2sum,coin_t);
 //----- Coin vs R-beta ----//
 hcoin_rbeta->Fill(rbeta,coin_t);
   if(cut_coin && cut_vz && cut_rbeta && cut_rpathl && cut_Rx){
    //---- Coin : AC -------//
     hcoin_ac1->Fill(Ra1sum,coin_t);
     hcoin_ac2->Fill(Ra2sum,coin_t);
   }
  }

 }// END FILL

 //=======================================================//
 //=============== AC Efficiency Study ===================//
 //=======================================================//
  int thmax=100;
  double ac1_th[thmax],ac2_th[thmax];
  double inte_ac1[thmax],inte_ac2[thmax],eff_ac1[100],eff_ac2[100];
  double bin_ac1max,bin_ac1min,bin_ac2max,bin_ac2min,bin_coinmin,bin_coinmax;
  for(int th=0;th<thmax;th++){
  ac1_th[th]=(max_ac1-min_ac1)/thmax*th+min_ac1;
  bin_ac1min=hcoin_ac1->GetXaxis()->FindBin(ac1_th[th]);
  bin_ac1max=hcoin_ac1->GetXaxis()->FindBin(max_ac1);
  bin_coinmin=hcoin_ac1->GetYaxis()->FindBin(coin_cutmin);
  bin_coinmax=hcoin_ac1->GetYaxis()->FindBin(coin_cutmax);
  inte_ac1[th]=hcoin_ac1->Integral(bin_ac1min,bin_ac1max,bin_coinmin,bin_coinmax);
  eff_ac1[th]=inte_ac1[th]/inte_ac1[0];
  //   cout<<"Efficiency of AC1 th:"<<th<<": "<<eff_ac1[th]<<endl;
  
  ac2_th[th]=(max_ac2-min_ac2)/thmax*th+min_ac2;
  bin_ac2min=hcoin_ac2->GetXaxis()->FindBin(ac2_th[th]);
  bin_ac2max=hcoin_ac2->GetXaxis()->FindBin(max_ac2);
  bin_coinmin=hcoin_ac2->GetYaxis()->FindBin(coin_cutmin);
  bin_coinmax=hcoin_ac2->GetYaxis()->FindBin(coin_cutmax);
  inte_ac2[th]=hcoin_ac2->Integral(bin_ac2min,bin_ac2max,bin_coinmin,bin_coinmax);
  eff_ac2[th]=inte_ac2[th]/inte_ac2[0];
  // cout<<"Efficiency of AC2 th:"<<th<<": "<<eff_ac2[th]<<endl;

  }

 //========= Draw ================//

TCanvas* ccoin=new TCanvas("ccoin","ccoin");
 ccoin->Divide(1,2);
 ccoin->cd(1);
 hcoin_t->Draw();
 ccoin->cd(2);
 hcoin_cut->Draw();

 TCanvas* cac=new TCanvas("cac","cac");
 cac->Divide(3,2);
 cac->cd(1);
 hbeta_ac1->Draw("colz");
 cac->cd(2);
 hbeta_ac2->Draw("colz");
 cac->cd(3);
 hac1_ac2->Draw("colz");
 cac->cd(4);
 hcoin_ac1->Draw("colz");
 cac->cd(5);
 hcoin_ac2->Draw("colz");
 cac->cd(6);
 hcoin_rbeta->Draw("colz");
 TCanvas* cpath=new TCanvas("cpath","cpath");
 cpath->Divide(2,2);
 cpath->cd(1);
 hcoin_rpathl->Draw("colz");
 cpath->cd(2);
 hcoin_lpathl->Draw("colz");
 cpath->cd(3);
 //hcoin_rpathlc->Draw("colz");
 cpath->cd(4);
 // hcoin_lpathlc->Draw("colz");

 
 TCanvas* ceff=new TCanvas("ceff","ceff");
 TGraph* ac1_eff=new TGraph(thmax,ac1_th,eff_ac1);
 TGraph* ac2_eff=new TGraph(thmax,ac2_th,eff_ac2);
 ac1_eff->SetTitle("Effieciency of AC1 ");
 ac1_eff->GetXaxis()->SetTitle("AC1 ADC Threshold [ch]");
 ac1_eff->GetYaxis()->SetTitle("Efficiency ");
 ac2_eff->SetTitle("Effieciency of AC2 ");
 ac2_eff->GetXaxis()->SetTitle("AC2 ADC Threshold [ch]");
 ac2_eff->GetYaxis()->SetTitle("Efficiency ");
 ceff->Divide(1,2);
 ceff->cd(1);
 ac1_eff->Draw();
 ceff->cd(2);
 ac2_eff->Draw();
 

 theApp->Run();
 return 0;
 } 
