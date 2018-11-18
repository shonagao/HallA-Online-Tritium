// Author K. Itabashi Aug. 30th
// HRS nnL experiment missing mass analysis

const double c=299792458;// [m/s]
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
   T->Add("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_111180.root");

 //============= Set Branch Status ==================//

  int max=10000; 
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
  double Rth[max],Rph[max],Lth[max],Lph[max];
  double Rbeta[max],Lbeta[max];
  double rs2pathl[max],rs0pathl[max],rtrpathl[max];
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
  T->SetBranchStatus("R.s2.trpath",rs2pathl); 
  T->SetBranchAddress("R.s2.trpath",rs2pathl); 
  T->SetBranchStatus("R.s0.trpath",rs0pathl); 
  T->SetBranchAddress("R.s0.trpath",rs0pathl); 
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
 T->SetBranchStatus("R.tr.beta",1);    
 T->SetBranchAddress("R.tr.beta",Rbeta); 


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

T->SetBranchStatus("L.tr.
beta",1);    
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

 //==================================================//
 // Time scale [ns] //
 // Energy scale [GeV] //

 double min_coin=-265;
 double max_coin=-240;
 int bin_coin=300;
 double min_beta=0.95;
 double max_beta=0.97;
 int bin_beta=6000;
 double min_adc=10.0;
 double max_adc=20000.;
 int bin_adc=max_adc-min_adc;
 TH1F* hmm=new TH1F("hmm","Missing mass Hist",5000,-0.1,2.);
 TH1F* hmm_c=new TH1F("hmm","Missing mass Hist",5000,-0.1,2.);
 TH1F* hmm1=new TH1F("hmm1","Missing mass Hist",5000,-0.1,2.);
 TH1F* hmm2=new TH1F("hmm2","Missing mass Hist",5000,-0.1,2.);
 TH1F* hmm3=new TH1F("hmm3","Missing mass Hist",5000,-0.1,2.);
 TH1F* hmm_ac=new TH1F("hmm_ac","Missing mass Hist with AC cut",500,-0.1,2.);
 TCanvas* c0=new TCanvas("c0","c0");
 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S0R-S0L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_t1=new TH1F("hcoin_t1","Coincidence time S0R-S0L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_t2=new TH1F("hcoin_t2","Coincidence time S0R-S0L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_t3=new TH1F("hcoin_t3","Coincidence time S0R-S0L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tbeta=new TH1F("hcoin_tbeta","Coincidence time S0R-S0L[sec] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_div=new TH1F("hcoin_div","Divede bin hcoin1/hcoin",bin_coin,min_coin,max_coin);
 TCanvas* c1=new TCanvas("c1","c1");
 TH2F* ha1_mm=new TH2F("ha1_mm","beta vs ac1 ADC sum hist",bin_beta,min_beta,max_beta,bin_adc,min_adc,max_adc);
 TCanvas* c3=new TCanvas("c3","c3");
 TH2F* ha2_mm=new TH2F("ha2_mm","beta vs ac2 ADC sum hist",bin_beta,min_beta,max_beta,bin_adc,min_adc,max_adc);
 TCanvas* cac=new TCanvas("cac","cac");
 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_adc,min_adc,max_adc,bin_adc,min_adc,max_adc);
 TCanvas* c4=new TCanvas("c4","c4");
 TH1F* hm2=new TH1F("hm2","Right ARM Mass Hist",500,-1.0,2.);
 TCanvas* c5=new TCanvas("c5","c5");

 double min_rpathl=0.0; double max_rpathl=30.; int bin_rpathl=10000;
 TH2F* hcoin_rpathl=new TH2F("hcoin_rpathl","Coinc time vs R Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
  
 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 double mh,mhc;
 double m2; 
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t;
 double rtof[16];
 double rbeta,rbeta_k;
 int i=8;
 int counts=0;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_trig,cut_beta;
 int nac1,nac2,n;
 double ac1_adc,ac2_adc;
 double tof_r,tof_l;
 double cos_ee_,cos_ek,cos_e_k,theta_r,theta_l,theta_rc,theta_lc;
 ac1_adc=100.;
 ac2_adc=2000.;
 double rpathl;
 TLine* lac1=new TLine(0,1.5,ac1_adc,ac1_adc);
 TLine* lac2=new TLine(0,1.5,ac2_adc,ac2_adc);

 for(int k=0;k<evnt;k++){
   T->GetEntry(k);
 
   cut_ac1=false;
   cut_ac2=false;
   cut_trig=true;
   cut_beta=false;
   if(Ra1sum<ac1_adc)cut_ac1=true;
   if(Ra2sum>ac2_adc)cut_ac2=true;
   if(rbeta>0.963 && rbeta<0.966)cut_beta=true;
   // if(trigger[0]==32)cut_trig=true;

   //  Lph[0]=Lph[0]-13.2*3.14/360;//rad
   // Rph[0]=Rph[0] +13.2*3.14/360;//rad
   // theta_r=atan(Rph[0])+13.2*3.14/360;
   // theta_l=atan(Lph[0])-13.2*3.14/360;
   // if(theta_l/(2*3.14)*360>12.0 && theta_l/(2*3.14)*360<12.4 && theta_r/(2*3.14)*360>-12.4 && theta_r/(2*3.14)*360<-12.0){
   /*  
 if(k<100){
     cout<<"theta_r w/o corr: "<<atan(Rph[0])/(2*3.14)*360.<<endl;
     cout<<"theta_l w/o corr: "<<atan(Lph[0])/(2*3.14)*360.<<endl;
     cout<<"theta_r w/ : "<<theta_r/(2*3.14)*360<<endl;
     cout<<"theta_l w/ : "<<theta_l/(2*3.14)*360<<endl;
     
     }*/


 pe_=Lp[0];//*sqrt(1+pow(Lth[0],2)+pow(Lph[0],2));
 pk =Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 ppi=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 
 theta_r=atan(Rph[0]);
 theta_l=atan(Lph[0]);   
 theta_rc=atan(Rph[0])-13.2*2*3.14/360;
 theta_lc=atan(Lph[0])+13.2*2*3.14/360;   
 cos_ee_=cos(theta_l);//1./sqrt(1+pow(Lph[0],2));
 cos_ek=cos(theta_r);//1./sqrt(1+pow(Rph[0],2));
 cos_e_k=cos(+theta_l-theta_r);//cos_ee_*cos_ek*(1-Lph[0]*Rph[0]);


 pe=hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Epi=sqrt(pow(ppi,2)+pow(mpi,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 mh=sqrt(pow(Ee+mtr-Ee_-Ek,2)-(pow(pe,2)+pow(pe_,2)+pow(pk,2)-2*pe*pe_*cos_ee_-2*pe*pk*cos_ek+2*pk*pe_*cos_e_k));
 mhc=sqrt(pow(Ee+mtr-Ee_-Ek,2)-(pow(pe,2)+pow(pe_,2)+pow(pk,2)-2*pe*pe_*cos(theta_lc)-2*pe*pk*cos(theta_rc)+2*pk*pe_*cos(theta_lc-theta_rc)));
 Ls2pads=Ls2tpads[0];
 Rs2pads=Rs2tpads[0];
 // tof_r=((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0-(-RF1[15]+RF1[9]))*tdc_time;
 // tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37])/2.0-(-RF1[47]+RF1[37]))*tdc_time;
 tof_r=((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0)*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37])/2.0)*tdc_time;
 coin_t=tof_r-tof_l;

 // coin_t= (-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0-(-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37])/2.0; //S2 coincidence time
 //coin_t= (RF1[43]-RF1[46]+RF1[44]-RF1[46])/2.0-(LF1[27]-LF1[30]-LF1[28]-LF1[30])/2.0; //S0 coincidence time
 rpathl=rtrpathl[0]-rs0pathl[0];


 // rtof[i]=(Rs2r_tc[i]+Rs2l_tc[i])/2.0-(Rs0r_tc[0]+Rs0l_tc[0])/2.0;
 // rbeta=(rs2pathl[0]-rs0pathl[0])/c/rtof[i];
 rbeta=pk/Ek;
 m2=sqrt((1./pow(Rbeta[0],2)-1)*pk);

 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2

 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0) && ( Ls0r_tc[0]>0 && Ls0l_tc[0]>0) && cut_trig){
hcoin_t->Fill(coin_t);
 hmm->Fill(mh);
 hmm_c->Fill(mhc);
 hcoin_rpathl->Fill(rpathl,coin_t);
 }
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0) && ( Ls0r_tc[0]>0 && Ls0l_tc[0]>0)&& cut_ac1 && cut_trig){
hcoin_t1->Fill(coin_t);
hmm1->Fill(mhc);
 }
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0) && ( Ls0r_tc[0]>0 && Ls0l_tc[0]>0)&& cut_ac2 && cut_trig){
hcoin_t2->Fill(coin_t);
hmm2->Fill(mhc);
 }
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0) && ( Ls0r_tc[0]>0 && Ls0l_tc[0]>0)&&cut_ac1 && cut_ac2 && cut_trig){
hcoin_t3->Fill(coin_t);
hmm3->Fill(mhc);
 }

 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0) && ( Ls0r_tc[0]>0 && Ls0l_tc[0]>0)&&cut_ac1 && cut_ac2 && cut_trig && cut_beta)hcoin_tbeta->Fill(coin_t);
 //if((Rs0r_tc>0 && Rs0l_tc>0) && ( Ls0r_tc>0 && Ls0l_tc>0))hmm->Fill(mh);
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0)&& (cut_ac1 && cut_ac2))hmm_ac->Fill(mh);
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0))ha1_mm->Fill(rbeta,Ra1sum);
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0))ha2_mm->Fill(rbeta,Ra2sum);
 if((Rs0r_tc[0]>0 && Rs0l_tc[0]>0))hm2->Fill(m2);


 }

   c0->cd();
 hmm->Draw();
 hmm_c->SetLineColor(1);
 hmm_c->Draw("same");
 /*
 hmm1->SetLineColor(2);
 hmm1->Draw("same");
 hmm2->SetLineColor(3);
 hmm2->Draw("same");
 hmm3->SetLineColor(1);
 hmm3->Draw("same");
 hmm_ac->SetLineColor(2);
 hmm_ac->Draw("same");
 */
 c1->cd();
 hcoin_t->Draw();
 hcoin_t1->SetLineColor(2);
 hcoin_t2->SetLineColor(3);
 hcoin_t3->SetLineColor(1);
 hcoin_tbeta->SetLineColor(4);
 hcoin_t1->Draw("same");
 hcoin_t2->Draw("same");
 hcoin_t3->Draw("same");
 hcoin_tbeta->Draw("same");

TCanvas *crpathl=new TCanvas("crpathl","crpathl");
crpathl->cd();
 hcoin_rpathl->Draw("colz");


 theApp->Run();
 return 0;
}

