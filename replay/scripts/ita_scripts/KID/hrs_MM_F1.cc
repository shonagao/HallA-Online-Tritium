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
const double mp=938.27e-3; //proton mass [GeV/c^2]
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
  int ch;
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_small.list";
  string ofname = "/pdf/test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:LR"))!=-1){
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
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      //      cout<<"-n : maximum number of analysed events"<<endl;
      //   cout<<"-p : input param file"<<endl;
      //     cout<<"-L : Left HRS"<<endl;
      //  cout<<"-R : Right HRS"<<endl;
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



  TApplication *theApp =new TApplication("App",&argc,argv);

  //============= Run Phase =====================//
  char* phase;
  phase="H_p1"; // Hydrogen run Phase1 111160-111220
 //phase="H_p2"; // Hydrogen run Phase1 111480-
 //phase="T_p1"; //Tritium run Phase1
 //=============================================//

 double tdc_time;
 if(phase=="H_p1")tdc_time=56.23e-3;//[ns]  Phase1 
 else if(phase=="H_p2")tdc_time=58e-3;//[ns] Phase2 
 else{cout<<"False Phase "<<endl; tdc_time=56.23e-3;}


//============ Opiton =============//



 if(draw_flag==0)gROOT->SetBatch(1);





 //========== Get Tree =============//


  
TChain*  T=new TChain("T");
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    // cout<<buf<<endl;
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

 // Target positon information //
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
  T->SetBranchStatus("L.tr.vz",1);    
  T->SetBranchAddress("L.tr.vz",Lvz);
 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double min_coin=-30;
 double max_coin=30.0;
 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 
 double min_coin_c=-15;
 double max_coin_c=5.0;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 

 double min_beta=0.95;
 double max_beta=0.97;
 int bin_beta=6000;

 double min_adc=-500.0;
 double max_adc=20000.;
 int bin_adc=max_adc-min_adc;

 double min_ac1=0.0;
 double max_ac1=5000.;
 int bin_ac1=max_ac1-min_ac1; 

 double min_ac2=0.0;
 double max_ac2=20000.;
 int bin_ac2=max_ac2-min_ac2; 

 //========= Coincidence Hist =======================//

 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t1=new TH1F("hcoin_t1","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t2=new TH1F("hcoin_t2","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t3=new TH1F("hcoin_t3","Coincidence time S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);

 //========= Path Correction Hist =================//

 double min_rpathl=28.0; double max_rpathl=30.; int bin_rpathl=500;
 double min_lpathl=28.0; double max_lpathl=30.; int bin_lpathl=500; 
 double min_coin_rc=80; double max_coin_rc=100;
 double bin_coin_rc=(max_coin_rc-min_coin_rc)/tdc_time;
        bin_coin_rc=(int)bin_coin_rc;
 double min_coin_lc=-120; double max_coin_lc=-100;
 double bin_coin_lc=(max_coin_lc-min_coin_lc)/tdc_time;
        bin_coin_lc=(int)bin_coin_lc;
 
 TH2F* hcoin_rpathl=new TH2F("hcoin_rpathl","Coinc time vs R Path Length Hist ",bin_rpathl,min_rpathl,max_rpathl,bin_coin,min_coin,max_coin); 
 hcoin_rpathl->SetTitle("R-Path Length vs Coinc-Time Hist w/o Path Correction; R-Path Length [m]; Coinc-time [ns]; Counts");
 TH2F* hcoin_lpathl=new TH2F("hcoin_lpathl","Coinc time vs L Path Length Hist ",bin_lpathl,min_lpathl,max_lpathl,bin_coin,min_coin,max_coin); 
hcoin_lpathl->SetTitle("L-Path Length vs Coinc-Time Hist w/o Path Correction; L-Path Length [m]; Coinc-time [ns]; Counts");
 TH2F* hcoin_rpathl_c=new TH2F("hcoin_rpathl_c","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_rpathl,max_rpathl,bin_coin_rc,min_coin_rc,max_coin_rc); 
 TH2F* hcoin_lpathl_c=new TH2F("hcoin_lpathl_c","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,bin_coin_lc,min_coin_lc,max_coin_lc);

 TH2F* hcoin_rpathl_cc=new TH2F("hcoin_rpathl_cc","Coinc time vs R Path Length Hist w/ correction",bin_rpathl,min_rpathl,max_rpathl,bin_coin_c,min_coin_c,max_coin_c);
hcoin_rpathl_cc->SetTitle("R-Path Length vs Coinc-Time Hist w/ Path Correction; R-Path Length [m]; Coinc-time [ns]; Counts");
 TH2F* hcoin_lpathl_cc=new TH2F("hcoin_lpathl_cc","Coinc time vs L Path Length Hist w/ correction",bin_lpathl,min_lpathl,max_lpathl,bin_coin_c,min_coin_c,max_coin_c);
  hcoin_lpathl_cc->SetTitle("L-Path Length vs Coinc-Time Hist w/ Path Correction; L-Path Length [m]; Coinc-time [ns]; Counts");

 //========= Coincidence time vs AC ==================//
 TH2F* hcoin_ac1=new TH2F("hcoin_ac1","Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 TH2F* hcoin_ac2=new TH2F("hcoin_ac2","Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);
 ha1_a2->SetTitle("ac1 vs ac2 ADC sum hist;AC1 adc;AC2 adc");
 
 //============= Missing Mass Hist ==================//
 double min_mm,max_mm;
 min_mm=-0.1; max_mm=0.2;
 int bin_mm=1000;
 TH1F* hmm=new TH1F("hmm","Missing mass Hist",bin_mm,min_mm,max_mm);
 TH2F* hmm_coin=new TH2F("hmm_coin","Missing mass vs Coin time  Hist",bin_coin_c,min_coin_c,max_coin_c,bin_mm,min_mm,max_mm);
 hmm_coin->SetTitle("Coincidence-time Cut in Missing Mass hist ;coinc-time [ns];Missing mass [GeV/c^2]; Counts");

 double min_z=-0.2; double max_z=0.2; int bin_z=1000;
 TH2F* hmm_vz=new TH2F("hmm_vz","Missing mass vs vs Z-vertex Hist",bin_z,min_z,max_z,bin_mm,min_mm,max_mm);
 hmm_vz->SetTitle("Z-vertex Cut in Missing mass ;Rz-vertex [m];Missing mass [GeV/c^2];Counts");


 int evnt=T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 cout<<"phase :"<<phase<<endl;
 cout<<"F1tdc time :"<<tdc_time<<endl;

 
 double mh,mtr,mm;
 if(phase=="H_p1" || phase=="H_p2"){
   mtr=mp;
   mh=ml;

   cout<<"Target mass: "<<mp<<endl;
   cout<<"hypernuclear mass: "<<mh<<endl;
  }


 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rtof[16];
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off; 
 double Rs2_tcorr,Ls2_tcorr;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double ac1_adc,ac2_adc;
 double tof_r,tof_l;
 double cos_ee_,cos_ek,cos_e_k,theta_r,theta_l,theta_rc,theta_lc; 
 
//===== AC Threshold =========//
  ac1_adc=150.;
  ac2_adc=1400.;
 //============================//
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr,pathl_off;

 double coin_offset;
 if(phase=="H_p1")coin_offset=-5.55;
 if(phase=="H_p2")coin_offset=20.85;
 //----- Cut Parameters ----------//
 double coin_cutmin=-1.0;
 double coin_cutmax=1.0; 
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

 //====== Missing Mass Parameters ===============//
 theta_r=atan(Rph[0]);

 theta_l=atan(Lph[0]);   
 theta_rc=atan(Rph[0])-13.2*2*3.14/360;
 theta_lc=atan(Lph[0])+13.2*2*3.14/360;   
 cos_ee_=cos(theta_l);//1./sqrt(1+pow(Lph[0],2));
 cos_ek=cos(theta_r);//1./sqrt(1+pow(Rph[0],2));
 cos_e_k=cos(+theta_l-theta_r);//cos_ee_*cos_ek*(1-Lph[0]*Rph[0]);

 

 mm=sqrt(pow(Ee+mtr-Ee_-Ek,2)-(pow(pe,2)+pow(pe_,2)+pow(pk,2)-2*pe*pe_*cos(theta_lc)-2*pe*pk*cos(theta_rc)+2*pk*pe_*cos(theta_lc-theta_rc)))-mh;

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
 coin_t=tof_r-tof_l-coin_offset; //coin time
 coin_tc=coin_t+rpath_corr-lpath_corr;//+pathl_off;


//====== Cut condition ========================// 
   cut_ac1=false;
   cut_ac2=false;
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_track=false;
   cut_s0=false;
   cut_coin=false;



   if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
   if(Ra1sum<ac1_adc)cut_ac1=true;
   if(Ra2sum>ac2_adc)cut_ac2=true;
   if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
   if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;

   if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
   if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
   if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
   if(coin_cutmin<coin_tc && coin_tc<coin_cutmax)cut_coin=true;
 //=======================================//




   //==========================================//
   //========= Fill Hist =====================//
   //========================================//


 if(cut_Rs2 && cut_Ls2 && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){

 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 hcoin_rpathl->Fill(rpathl,coin_t);
 hcoin_rpathl_c->Fill(rpathl,coin_t+rpath_corr);
 hcoin_rpathl_cc->Fill(rpathl,coin_tc);
 hcoin_lpathl->Fill(lpathl,coin_t);
 hcoin_lpathl_c->Fill(lpathl,coin_t-lpath_corr);
 hcoin_lpathl_cc->Fill(lpathl,coin_tc);


 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2
 }
 if(cut_Rs2 && cut_Ls2 && cut_ac1 && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t1->Fill(coin_tc); // AC1 cut
   hcoin_ac2->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
 }
 if(cut_Rs2 && cut_Ls2 && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t2->Fill(coin_tc); //AC2 cut
   hcoin_ac1->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut) 
 }

 if(cut_Rs2 && cut_Ls2  && cut_ac1 && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track && cut_s0){
   hcoin_t3->Fill(coin_tc); //AC1 & AC2 cut
   hmm_coin->Fill(coin_tc,mm); // Missing mass
   if(cut_coin)hmm->Fill(mm); // Coin vs Misiing mass
 
 }


 if(cut_Rs2 && cut_Ls2  && cut_ac1 && cut_ac2 && cut_lpathl && cut_rpathl && cut_track && cut_s0)hmm_vz->Fill(Rvz[0],mm); // Coin vs Z-vertex


 }



 n=hcoin_t->GetEntries();
 nac1=hcoin_t1->GetEntries();
 nac2=hcoin_t2->GetEntries();
 nac3=hcoin_t3->GetEntries();
  


 //============== Efficiency analysis =======================//

 
 int iter_ac1=50; int iter_ac2=100; int iter_max=iter_ac1+iter_ac2;

 TF1* facc[iter_max][2];
 TF1* fpi[iter_max][2];
 TF1* fk[iter_max][2];
 TF1* fcoin[iter_max][2];
 TF1* fp[iter_max][2];
 TH1D* hcoin_ac1_p[iter_max];
 TH1D* hcoin_ac2_p[iter_max]; 
 TF1* facc_t1_def;
 TF1* fpi_t1_def;

 TF1* fk_t1_def;
 TF1* fp_t1_def;
 TF1* facc_t2_def;

 TF1* fpi_t2_def;
 TF1* fk_t2_def;
 TF1* fp_t2_def;
 TF1* facc_t3_def;
 TF1* fpi_t3_def;
 TF1* fk_t3_def;
 TF1* fp_t3_def;
 TGraph* gsum_pi_ac1=new TGraph();
 TGraph* gsum_p_ac1=new TGraph();
 TGraph* gsum_k_ac1=new TGraph();
 TGraph* grate_k_ac1=new TGraph();
 TGraph* grate_p_ac1=new TGraph();
 TGraph* grate_pi_ac1=new TGraph();
 TGraph* gsum_pi_ac2=new TGraph();
 TGraph* gsum_p_ac2=new TGraph();

 TGraph* gsum_k_ac2=new TGraph();
 TGraph* grate_k_ac2=new TGraph();


 TGraph* grate_pi_ac2=new TGraph();
 TGraph* grate_p_ac2=new TGraph();
 //--- TGraph Setting ----//
 gsum_pi_ac1->SetTitle("SUM of Pion vs AC1 Threshold;AC1 ADC-th [ch];Pion Events [Counts]");
 gsum_p_ac1->SetTitle("SUM of Proton vs AC1 Threshold;AC1 ADC-th [ch];Proton Events [Counts]");
 gsum_k_ac1->SetTitle("SUM of Kaon vs AC1 Threshold;AC1 ADC-th [ch];Kaon Events [Counts]");
 grate_k_ac1->SetTitle("Kaoin Survival rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Survival rate ");
 grate_pi_ac1->SetTitle("Pion Survuval rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Surival rate ");
 grate_p_ac1->SetTitle("Proton Survival rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Survival rate ");
 gsum_pi_ac1->SetMarkerStyle(21);
 gsum_pi_ac1->SetMarkerColor(kRed);
 gsum_pi_ac1->SetMarkerSize(0.5);
 gsum_p_ac1->SetMarkerStyle(21);
 gsum_p_ac1->SetMarkerColor(kRed);
 gsum_p_ac1->SetMarkerSize(0.5);
 gsum_k_ac1->SetMarkerStyle(21);
 gsum_k_ac1->SetMarkerColor(kRed);
 gsum_k_ac1->SetMarkerSize(0.5);
 grate_k_ac1->SetMarkerStyle(21);
 grate_k_ac1->SetMarkerColor(kBlue);
 grate_k_ac1->SetMarkerSize(0.5);
 grate_p_ac1->SetMarkerStyle(21);
 grate_p_ac1->SetMarkerColor(kBlue);
 grate_p_ac1->SetMarkerSize(0.5);
 grate_pi_ac1->SetMarkerStyle(21);
 grate_pi_ac1->SetMarkerColor(kBlue);

 grate_pi_ac1->SetMarkerSize(0.5);

 gsum_pi_ac2->SetTitle("SUM of Pion vs AC2 Threshold;AC2 ADC-th [ch];Pion Events [Counts]");
 gsum_p_ac2->SetTitle("SUM of Proton vs AC2 Threshold;AC2 ADC-th [ch];Proton Events [Counts]");
 gsum_k_ac2->SetTitle("SUM of Kaon vs AC2 Threshold;AC2 ADC-th [ch];Kaon Events [Counts]");
 grate_k_ac2->SetTitle("Kaoin S/N rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_pi_ac2->SetTitle("Pion Survuval rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_p_ac2->SetTitle("Proton Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 gsum_pi_ac2->SetMarkerStyle(21);
 gsum_pi_ac2->SetMarkerColor(kRed);
 gsum_pi_ac2->SetMarkerSize(0.5);
 gsum_p_ac2->SetMarkerStyle(21);
 gsum_p_ac2->SetMarkerColor(kRed);
 gsum_p_ac2->SetMarkerSize(0.5);
 gsum_k_ac2->SetMarkerStyle(21);
 gsum_k_ac2->SetMarkerColor(kRed);
 gsum_k_ac2->SetMarkerSize(0.5);
 grate_k_ac2->SetMarkerStyle(21);
 grate_k_ac2->SetMarkerColor(kBlue);
 grate_k_ac2->SetMarkerSize(0.5);
 grate_p_ac2->SetMarkerStyle(21);
 grate_p_ac2->SetMarkerColor(kBlue);
 grate_p_ac2->SetMarkerSize(0.5);
 grate_pi_ac2->SetMarkerStyle(21);
 grate_pi_ac2->SetMarkerColor(kBlue);
 grate_pi_ac2->SetMarkerSize(0.5);
 
//--- Parameters -----//

 double kmin[iter_max][2],kmax[iter_max][2];
 double inte_ktot[iter_max][2], inte_ksig[iter_max][2];
 double p0_acc[iter_max][2], p1_acc[iter_max][2];
 double n_p[iter_max][2],sig_p[iter_max][2],mean_p[iter_max][2];
 double n_pi[iter_max][2],sig_pi[iter_max][2],mean_pi[iter_max][2];

 double n_k[iter_max][2],sig_k[iter_max][2],mean_k[iter_max][2];
 int bin_ac1_adc,bin_min_ac1,bin_ac2_adc,bin_max_ac2,bin_min_ac2;
 double sum_k[iter_max][2],sum_p[iter_max][2],sum_pi[iter_max][2]; 
 double inte_acc[iter_max][2];
 double th_ac1,th_ac2;
 int bin_th_ac1,bin_th_ac2; 
 double nk[iter_max][2],npi[iter_max][2],np[iter_max][2];
 double max_nk[2],max_npi[2],max_np[2];
//---- Defolt parameters -----------//
 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_mean_k,def_sig_k,def_acc;

//==== HYdrogen Run Phase1 111160-111220========//



 if(phase=="H_p1"){
 def_sig_p=0.852;  def_mean_p=-13.82-coin_offset;
 def_sig_pi=0.443; def_mean_pi=-2.63-coin_offset;
 def_sig_k=0.644;  def_mean_k=-5.55-coin_offset;
 def_acc=27.7; 
 }else if(phase=="H_p2"){


 }



 double def_t2_k,def_t2_p,def_t2_pi,def_t1_k,def_t1_pi,def_t1_p,def_t3_k,def_t3_pi,def_t3_p;
 double t1sum_k,t1sum_pi,t1sum_p,t2sum_k,t2sum_pi,t2sum_p,t3sum_p,t3sum_pi,t3sum_k;



 int ac1_p,ac1_pi,ac1_k,ac2_k,ac2_p,ac2_pi;
 ac1_p=0;ac1_pi=0;ac1_k=0;ac2_p=0;ac2_pi=0;ac2_k=0;
 double t1sig_k,t2sig_k,t1mean_k,t2mean_k,t3sig_k,t3mean_k;
 double t1sig_p,t2sig_p,t1mean_p,t2mean_p,t3sig_p,t3mean_p;
 double t1sig_pi,t2sig_pi,t1mean_pi,t2mean_pi,t3sig_pi,t3mean_pi;
 double min_k,max_k,min_p,max_p,min_pi,max_pi;
 min_k=def_mean_k-2*def_sig_k;  max_k=def_mean_k+2*def_sig_k;
 min_p=def_mean_p-2*def_sig_p;  max_p=def_mean_p+2*def_sig_p;
 min_pi=def_mean_pi-2*def_sig_pi;  max_pi=def_mean_pi+2*def_sig_pi;

 //=============== Maximum Event Number Analysis =======================================//
 
 facc_t1_def=new TF1("facc_t1_def","[0]",min_coin_c,max_coin_c);

 fk_t1_def=new TF1("fk_t1_def","gausn",min_coin_c,max_coin_c);
 fpi_t1_def=new TF1("fpi_t1_def","gausn",min_coin_c,max_coin_c);
 fp_t1_def=new TF1("fp_t1_def","gausn",min_coin_c,max_coin_c);
 facc_t2_def=new TF1("facc_t2_def","[0]",min_coin_c,max_coin_c);

 fk_t2_def=new TF1("fk_t2_def","gausn",min_coin_c,max_coin_c);
 fpi_t2_def=new TF1("fpi_t2_def","gausn",min_coin_c,max_coin_c);
 fp_t2_def=new TF1("fp_t2_def","gausn",min_coin_c,max_coin_c);
 facc_t3_def=new TF1("facc_t3_def","[0]",min_coin_c,max_coin_c);
 fk_t3_def=new TF1("fk_t3_def","gausn",min_coin_c,max_coin_c);
 fpi_t3_def=new TF1("fpi_t3_def","gausn",min_coin_c,max_coin_c);
 fp_t3_def=new TF1("fp_t3_def","gausn",min_coin_c,max_coin_c);

 TF1* fcoin_t1 =new TF1("fcoin_t1","gausn(0)+gausn(3)+gausn(6)+pol1(9)",min_coin_c,max_coin_c);

 fcoin_t1->SetParameters(1115,def_mean_p,def_sig_p,140.,def_mean_pi,def_sig_pi,49.7,def_mean_k,def_sig_k,def_acc); 
 TF1* fcoin_t2 =new TF1("fcoin_t2","gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin_t2->SetParameters(1115,def_mean_p,def_sig_p,140.,def_mean_pi,def_sig_pi,49.7,def_mean_k,def_sig_k,def_acc); 
 TF1* fcoin_t3 =new TF1("fcoin_t3","gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin_t3->SetParameters(1115,def_mean_p,def_sig_p,140.,def_mean_pi,def_sig_pi,49.7,def_mean_k,def_sig_k,def_acc);  

 
 fcoin_t1->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t1->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
 fcoin_t1->SetParLimits(4,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t1->SetParLimits(5,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t1->SetParLimits(7,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t1->SetParLimits(8,0.8*def_sig_k,1.2*def_sig_k);
 
 fcoin_t2->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t2->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
 fcoin_t2->SetParLimits(4,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t2->SetParLimits(5,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t2->SetParLimits(7,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t2->SetParLimits(8,0.8*def_sig_k,1.2*def_sig_k);

 fcoin_t3->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t3->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
 fcoin_t3->SetParLimits(4,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t3->SetParLimits(5,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t3->SetParLimits(7,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t3->SetParLimits(8,0.8*def_sig_k,1.2*def_sig_k);
 



  hcoin_t1->Fit("fcoin_t1","R","",min_coin_c,max_coin_c);
  hcoin_t2->Fit("fcoin_t2","R","",min_coin_c,max_coin_c);
  hcoin_t3->Fit("fcoin_t3","R","",min_coin_c,max_coin_c);
 
 def_t1_p=fcoin_t1->GetParameter(0); t1mean_p=fcoin_t1->GetParameter(1); t1sig_p=fcoin_t1->GetParameter(2); 
 def_t1_pi=fcoin_t1->GetParameter(3);t1mean_pi=fcoin_t1->GetParameter(4); t1sig_pi=fcoin_t1->GetParameter(5);  
 def_t1_k=fcoin_t1->GetParameter(6);t1mean_k=fcoin_t1->GetParameter(7); t1sig_k=fcoin_t1->GetParameter(8);

 def_t2_p=fcoin_t2->GetParameter(0); t2mean_p=fcoin_t2->GetParameter(1); t2sig_p=fcoin_t2->GetParameter(2); 
 def_t2_pi=fcoin_t2->GetParameter(3);t2mean_pi=fcoin_t2->GetParameter(4); t2sig_pi=fcoin_t2->GetParameter(5);  
 def_t2_k=fcoin_t2->GetParameter(6);t2mean_k=fcoin_t2->GetParameter(7); t2sig_k=fcoin_t2->GetParameter(8); 


 

 def_t3_p=fcoin_t3->GetParameter(0); t3mean_p=fcoin_t3->GetParameter(1); t3sig_p=fcoin_t3->GetParameter(2); 
 def_t3_pi=fcoin_t3->GetParameter(3);t3mean_pi=fcoin_t3->GetParameter(4); t3sig_pi=fcoin_t3->GetParameter(5);  
 def_t3_k=fcoin_t3->GetParameter(6);t3mean_k=fcoin_t3->GetParameter(7); t3sig_k=fcoin_t3->GetParameter(8); 


 
 fk_t1_def->SetParameters(def_t1_k,t1mean_k,t1sig_k);
 fpi_t1_def->SetParameters(def_t1_pi,t1mean_pi,t1sig_pi);
 fp_t1_def->SetParameters(def_t1_p,t1mean_p,t1sig_p);
 fk_t2_def->SetParameters(def_t2_k,t2mean_k,t2sig_k);
 fpi_t2_def->SetParameters(def_t2_pi,t2mean_pi,t2sig_pi);
 fp_t2_def->SetParameters(def_t2_p,t2mean_p,t2sig_p);
 

 fk_t3_def->SetParameters(def_t3_k,t3mean_k,t3sig_k);
 fpi_t3_def->SetParameters(def_t3_pi,t3mean_pi,t3sig_pi);
 fp_t3_def->SetParameters(def_t3_p,t3mean_p,t3sig_p);


 t1sum_k=def_t1_k/tdc_time;   t2sum_k=def_t2_k/tdc_time;  
 t1sum_p=def_t1_p/tdc_time;   t2sum_p=def_t2_p/tdc_time;   
 t1sum_pi=def_t1_pi/tdc_time; t2sum_pi=def_t2_pi/tdc_time; 
 t3sum_k=def_t3_k/tdc_time;
 t3sum_p=def_t3_p/tdc_time;

 t3sum_pi=def_t3_pi/tdc_time;

 max_nk[0]=0.0; max_nk[1]=0.0;
 max_np[0]=0.0; max_np[1]=0.0;

 max_npi[0]=0.0; max_npi[1]=0.0;






//----------------------------------//
 



 //=======================================================================//
 //==================== Threshold Change =================================//
 //=======================================================================//

 /*

     for(int i=0;i<iter_ac1;i++){

       th_ac1=min_ac1+(ac1_adc-min_ac1)/iter_ac1*i;
       th_ac2=min_ac2+(ac2_adc-min_ac2)/iter_ac1*i;	
 bin_th_ac1=hcoin_ac1->GetYaxis()->FindBin(th_ac1);
 bin_min_ac1=hcoin_ac1->GetYaxis()->FindBin(min_ac1);
 hcoin_ac1_p[i]=hcoin_ac1->ProjectionX(Form("hcoin_ac1_p[%d]",i),bin_min_ac1,bin_th_ac1);
 
 bin_th_ac2=hcoin_ac2->GetYaxis()->FindBin(th_ac2);

 bin_min_ac2=hcoin_ac2->GetYaxis()->FindBin(min_ac2);
 hcoin_ac2_p[i]=hcoin_ac2->ProjectionX(Form("hcoin_ac2_p[%d]",i),bin_min_ac2,bin_th_ac2);




 / /--- Initial Parameters -----------//


 facc[i][1]=new TF1(Form("facc[%d][1]",i),"[0]",min_coin_c,max_coin_c); 
 facc[i][0]=new TF1(Form("facc[%d][0]",i),"[0]",min_coin_c,max_coin_c); 
 hcoin_ac1_p[i]->Fit(Form("facc[%d][0]",i),"R","",min_coin_c,min_coin_c+2.5);

 p0_acc[i][0]=facc[i][0]->GetParameter(0);
 hcoin_ac2_p[i]->Fit(Form("facc[%d][1]",i),"R","",min_coin_c,min_coin_c+2.5);
 p0_acc[i][1]=facc[i][1]->GetParameter(0);


 //------- AC1 --------------// 


 fp[i][0]=new TF1(Form("fp[%d][0]",i),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi[i][0] =new TF1(Form("fpi[%d][0]",i),"gausn+pol0(3)",min_coin_c,max_coin_c);
 fk[i][0]=new TF1(Form("fk[%d][0]",i),"gausn+pol0(3)",min_coin_c,max_coin_c);
  fp[i][0]->SetParameter(1,def_mean_p); fp[i][0]->SetParameter(2,def_sig_p); fp[i][0]->SetParameter(3
,p0_acc[i][0]);
 fpi[i][0]->SetParameter(1,def_mean_pi); fpi[i][0]->SetParameter(2,def_sig_pi); fpi[i][0]->SetParameter(3,p0_acc[i][0]);

 
 fgk[i][0]->SetParameter(1,def_mean_k); fk[i][0]->SetParameter(2,def_sig_k); fk[i][0]->SetParameter(3,p0_acc[i][0]);





 fcoin[i][0] =new TF1(Form("fcoin[%d][0]",i),"gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin[i][0]->SetTitle(Form("Coin-Time w AC cut  (AC1<%d ch && AC2>%d ch);Coin time [ns];Counts [1/56 ns]",th_ac1,ac2_adc));
 //------- AC2 -------------//

 fp[i][1]=new TF1(Form("fp[%d][1]",i),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi[i][1] =new TF1(Form("fpi[%d][1]",i),"gausn+pol0(3)",min_coin_c,max_coin_c);

 fk[i][1]=new TF1(Form("fk[%d][1]",i),"gausn+pol0(3)",min_coin_c,max_coin_c);

 fp[i][1]->SetParameter(1,def_mean_p); fp[i][1]->SetParameter(2,def_sig_p); fp[i][1]->SetParameter(3,p0_acc[i][1]);
 fpi[i][1]->SetParameter(1,def_mean_pi); fpi[i][1]->SetParameter(2,def_sig_pi); fpi[i][1]->SetParameter(3,p0_acc[i][1]);
 fk[i][1]->SetParameter(1,def_mean_k); fk[i][1]->SetParameter(2,def_sig_k); fk[i][1]->SetParameter(3,

p0_acc[i][1]);
 fcoin[i][1] =new TF1(Form("fcoin[%d][1]",i),"gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin[i][1]->SetTitle(Form("Coin-Time w AC cut (AC1 Cut<%d ch && AC2 Cut>%d ch


);Coin time [ns];Counts [1/56 ns]",ac1_adc,th_ac2));
 //----------------------------------//
 
 //----- AC1 Fitting -----------// 
 hcoin_ac1_p[i]->Fit(Form("fp[%d][0]",i),"R","",t1mean_p-3*t1sig_p,t1mean_p+3*t1sig_p);
 n_p[i][0]=fp[i][0]->GetParameter(0);
 mean_p[i][0]=fp[i][0]->GetParameter(1);

 sig_p[i][0]=fp[i][0]->GetParameter(2);
 hcoin_ac1_p[i]->Fit(Form("fpi[%d][0]",i),"R","",t1mean_pi-3*t1sig_pi,t1mean_pi+3*t1sig_pi);

 n_pi[i][0]=fpi[i][0]->GetParameter(0);
 mean_pi[i][0]=fpi[i][0]->GetParameter(1);
 sig_pi[i][0]=fpi[i][0]->GetParameter(2);
 hcoin_ac1_p[i]->Fit(Form("fk[%d][0]",i),"R","",t1mean_k-3*t1sig_k,t1mean_k+3*t1sig_k);
 n_k[i][0]=fk[i][0]->GetParameter(0);
 mean_k[i][0]=fk [i][0]->GetParameter(1);
 sig_k[i][0]=fk[i][0]->GetParameter(2);


  //----- AC2 Fitting -----------// 
 hcoin_ac2_p[i]->Fit(Form("fp[%d][1]",i),"R","",t2mean_p-3*t2sig_p,t2mean_p+3*t2sig_p);
 n_p[i][1]=fp[i][1]->GetParameter(0);
 mean_p[i][1]=fp[i][1]->GetParameter(1);
 sig_p[i][1]=fp[i][1]->GetParameter(2);



 hcoin_ac2_p[i]->Fit(Form("fpi[%d][1]",i),"R","",t2mean_pi-3*t2sig_pi,t2mean_pi+3*t2sig_pi);
 n_pi[i][1]=fpi[i][1]->GetParameter(0);
 mean_pi[i][1]=fpi[i][1]->GetParameter(1);
 sig_pi[i][1]=fpi[i][1]->GetParameter(2);
 hcoin_ac2_p[i]->Fit(Form("fk[%d][1]",i),"R","",t2mean_k-3*t2sig_k,t2mean_k+3*t2sig_k);
 n_k[i][1]=fk[i][1]->GetParameter(0);
 mean_k[i][1]=fk[i][1]->GetParameter(1);
 sig_k[i][1]=fk[i][1]->GetParameter(2);

 //----- AC1 Coint Fitting ---------//


 fcoin[i][0]->SetParameters(n_pi[i][0],mean_pi[i][0],sig_pi[i][0],n_k[i][0],mean_k[i][0],sig_k[i][0],n_p[i][0],mean_p[i][0],sig_p[i][0],p0_acc[i][0]);
 fcoin[i][0]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin[i][0]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);

 fcoin[i][0]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);

 fcoin[i][0]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin[i][0]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin[i][0]->S etParLimits(8,0.8*def_sig_p,1.2*def_sig_p);
 fcoin[i][0]->SetParLimits(9,0.8*p0_acc[i][0],1.2*p0_acc[i][0]);
 hcoin_ac1_p[i]->Fit(Form("fcoin[%d][0]",i),"R","",min_coin_c,max_coin_c);
 n_pi[i][0]=fcoin[i][0]->GetParameter(0); mean_pi[i][0]=fcoin[i][0]->GetParameter(1);  sig_pi[i][0]=fcoin[i][0]->GetParameter(2);

 n_k[i][0]=fcoin[i][0]->GetParameter(3); mean_k[i][0]=fcoin[i][0]->GetParameter(4);  sig_k[i][0]=fcoin[i][0]->GetParameter(5);
 n_p[i][0]=fcoin[i][0]->GetParameter(6); mean_p[i][0]=fcoin[i][0]->GetParameter(7);  sig_p[i][0]=fcoin[i][0]->GetParameter(8);
 p0_acc[i][0]=fcoin[i][0]->GetParameter(9);



fp[i][0]->SetParameter(0,n_p[i][0]); fp[i][0]->SetParameter(1,mean_p[i][0]);fp[i][0]->SetParameter(2,sig_p[i][0]); fp[i][0]->SetParameter(3,p0_acc[i][0]);

fpi[i][0]->SetParameter(0,n_pi[i][0]); fpi[i][0]->SetParameter(1,mean_pi[i][0]);fpi[i][0]->SetParameter(2,sig_pi[i][0]); fpi[i][0]->SetParameter(3,p0_acc[i][0]);
fk[i][0]->SetParameter(0,n_k[i][0]); fk[i][0]->SetParameter(1,mean_k[i][0]);fk[i][0]->SetParameter(2,sig_k[i][0]); fk[i][0]->SetParameter(3,p0_acc[i][0]);


 //----- AC2 Coint Fitting ---------//

 fcoin[i][1]->SetParameters(n_pi[i][1],mean_pi[i][1],sig_pi[i][1],n_k[i][1],mean_k[i][1],sig_k[i][1],n_p[i][1],mean_p[i][1],sig_p[i][1],p0_acc[i][1]);
 fcoin[i][1]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin[i][1]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin[i][1]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin[i][1]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin[i][1]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin[i][1]->SetParLimits(8,0.8*def_sig_p,1.2*def_sig_p);
 fcoin[i][1]->SetParLimits(9,0.8*p0_acc[i][1],1.2*p0_acc[i][1]);
 hcoin_ac2_p[i]->Fit(Form("fcoin[%d][1]",i),"R","",min_coin_c,max_coin_c);

 n_pi[i][1]=fcoin[i][1]->GetParameter(0); mean_pi[i][1]=fcoin[i][1]->GetParameter(1);  sig_pi[i][1]=fcoin[i][1]->GetParameter(2);
 n_k[i][1]=fcoin[i][1]->GetParameter(3); mean_k[i][1]=fcoin[i][1]->GetParameter(4);  sig_k[i][1]=fcoin[i][1]->GetParameter(5);
 n_p[i][1]=fcoin[i][1]->GetParameter(6); mean_p[i][1]=fcoin[i][1]->GetParameter(7);  sig_p[i][1]=fcoi
n[i][1]->GetParameter(8);
 p0_acc[i][1]= fcoin[i][1]->GetParameter(9);

fp[i][1]->SetParameter(0,n_p[i][1]); fp[i][1]->SetParameter(1,mean_p[i][1]);fp[i][1]->SetParameter(2,sig_p[i][1]); fp[i][1]->SetParameter(3,p0_acc[i][1]);
fpi[i][1]->SetParameter(0,n_pi[i][1]); fpi[i][1]->SetParameter(1,mean_pi[i][1]);fpi[i][1]->SetParameter(2,sig_pi[i][1]); fpi[i][1]->SetParameter(3,p0_acc[i][1]);
fk[i][1]->SetParameter(0,n_k[i][1]); fk[i][1]->SetParameter(1,mean_k[i][1]);fk[i][1]->SetParameter(2,sig_k[i][1]); fk[i][1]->SetParameter(3,p0_acc[i][1]);

 
  //---- AC1 ----//
 sum_k[i][0]=n_k[i][0]/tdc_time; 
 sum_pi[i][0]=n_pi[i][0]/tdc_time;
 sum_p[i][0]=n_p[i][0]/tdc_time;



 //---- AC2 ----//
 sum_k[i][1]=n_k[i][1]/tdc_time; 
 sum_pi[i][1] =n_pi[i][1]/tdc_time;
 sum_p[i][1]=n_p[i][1]/tdc_time;
 //---- AC1 ----//
 gsum_pi_ac1->SetPoint(i,th_ac1,sum_pi[i][0]);

 gsum_p_ac1->SetPoint(i,th_ac1,sum_p[i][0]);
 gsum_k_ac1->SetPoint(i,th_ac1,sum_k[i][0]);
 grate_k_ac1->SetPoint(i,th_ac1,sum_k[i][0]/t3sum_k);
 grate_p_ac1->SetPoint(i,th_ac1,sum_p[i][0]/t3sum_p);

 grate_pi_ac1->SetPoint(i,th_ac1,sum_pi[i][0]/t3sum_pi);
 //---- AC2 ----//

 gsum_pi_ac2->SetPoint(i,th_ac2,sum_pi[i][1]);
 gsum_p_ac2->SetPoint(i,th_ac2,sum_p[i][1]);
 gsum_k_ac2->SetPoint(i,th_ac2,sum_k[i][1]);
 grate_k_ac2->SetPoint(i,th_ac2,sum_k[i][1]/t3sum_k);
 grate_p_ac2->SetPoint(i,th_ac2,sum_p[i][1]/t3sum_p);
 grate_pi_ac2->SetPoint(i,th_ac2,sum_pi[i][1]/t3sum_pi);

 
 

 if(n_k[i][0 ]>max_nk[0]){max_nk[0]=n_k[i][0]; ac1_k=i;}   

 if(n_k[i][1]>max_nk[1]){max_nk[1]=n_k[i][1]; ac2_k=i;}   
 if(n_p[i][0]>max_np[0]){max_np[0]=n_p[i][0]; ac1_p=i;}   

 if(n_p[i][1]>max_np[1]){max_np[1]=n_p[i][1]; ac2_p=i;}
 if(n_pi[i][0]>max_npi[0]){max_npi[0]=n_pi[i][0]; ac1_pi=i;}
 mif(n_pi[i][1]>max_npi[1]){max_npi[1]=n_pi[i][1]; ac2_pi=i;}
     

}

  

 */



 //======== Draw TCanvas ==============//
  



 
 fk_t1_def->SetLineColor(4);
 fk_t1_def->SetFillColor(4); 
 fk_t1_def->SetFillStyle(3001);
 fpi_t1_def->SetLineColor(46); 
 fpi_t1_def->SetFillStyle(3001);
 fpi_t1_def->SetFillColor(46);
 fp_t1_def->SetFillStyle(3001);
 fp_t1_def->SetLineColor(8);
 fp_t1_def->SetFillColor(8);
 fk_t2_def->SetFillStyle(3001);
 fk_t2_def->SetLineColor(4);
 fk_t2_def->SetFillColor(4);
 fpi_t2_def->SetFillStyle(3001);
 fpi_t2_def->SetLineColor(46);
 fpi_t2_def->SetFillColor(46);
 fp_t2_def->SetFillStyle(3001);
 fp_t2_def->SetLineColor(8);
 fp_t2_def->SetFillColor(8);  
 
 fk_t3_def->SetFillStyle(3001);
 fk_t3_def->SetLineColor(4);
 fk_t3_def->SetFillColor(4);
 fpi_t3_def->SetFillStyle(3001);
 fpi_t3_def->SetLineColor(46);
 fpi_t3_def->SetFillColor(46);
 fp_t3_def->SetFillStyle(3001);
 fp_t3_def->SetLineColor(8);
 fp_t3_def->SetFillColor(8);  
 

 TLine* lac=new TLine(min_coin_c,ac1_adc,max_coin_c,ac1_adc);
 lac->SetLineWidth(2);
 lac->SetLineColor(2);

 TCanvas* cmm_cut=new TCanvas("cmm_cut","cmm_cut");
 cmm_cut->Divide(1,2);
 cmm_cut->cd(1);
 hmm_coin->Draw("colz");
 lac->DrawLine(coin_cutmin,min_mm,coin_cutmin,max_mm); 
 lac->DrawLine(coin_cutmax,min_mm,coin_cutmax,max_mm); 
 cmm_cut->cd(2);
 hmm_vz->Draw("colz"); 
 lac->DrawLine(Rvz_cutmin,min_mm, Rvz_cutmin,max_mm); 
 lac->DrawLine(Rvz_cutmax,min_mm, Rvz_cutmax,max_mm); 

 

TCanvas* ccoin=new TCanvas("ccoin","ccoin");
 ccoin->cd(1);
 hcoin_t3->Draw();
 fcoin_t3->Draw("same");
 fk_t3_def->Draw("same");
 fp_t3_def->Draw("same");
 fpi_t3_def->Draw("same");

 //ccoin->Divide(1,3);
 /* 
ccoin->cd(1);
 hcoin_t1->Draw();  
 fcoin_t1->Draw("same");
 fk_t1_def->Draw("same");
 fp_t1_def->Draw("same");
 fpi_t1_def->Draw("same");
 ccoin->cd(2);
 hcoin_t2->Draw();
 hcoin_t2->Draw();  
 fcoin_t2->Draw("same");
 fk_t2_def->Draw("same");
 fp_t2_def->Draw("same");
 fpi_t2_def->Draw("same");
 */ 

 /*
 TCanvas* c _pathc =new TCanvas("c_pathc","c_pathc");
 c_pathc->Divide(1,2);
 c_pathc->cd(1);
 hcoin_t->Draw();
 c_pathc-> cd(2);
 hcoin_tc->Draw();
 */

 TCanvas *crpathl=new TCanvas("crpathl","crpathl");
 crpathl->Divide(2,2);
 crpathl->cd(1);
 hcoin_rpathl->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);

 crpathl->cd(2);
 hcoin_lpathl->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 
crpathl->cd(3);
 hcoin_rpathl_cc->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);
 crpathl->cd(4);
 hcoin_lpathl_cc ->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 


 /*
 TCanvas* ccoin=new TCanvas("ccoin","ccoin");
 ccoin->Divide(1,3);

 ccoin->cd(1);
 hcoin_t1->Draw();  
 fp_t1_def->Draw("same");
 fk_t1_def->Draw("same");
 fpi_t1_def->Draw("same");
 fcoin_t1->Draw("same");
 ccmoin->cd(2);
 hcoin_t2->Draw();
 fcoin_t2->Draw("same");
 fk_t2_def->Draw("same");

 ccoin->cd(3);
 hcoin_t3->Draw();
 fp_t3_def->Draw("same");

 fpi_t3_def->Draw("same");
 fcoin_t3->Draw("same"); 

 */
 
 TCanvas* ccoin_ac=new TCanvas("ccoin_ac","ccoin_ac");

 ccoin_ac->Divide(1,2);
 ccoin_ac->cd(1);
 hcoin_ac1->Draw("colz");
 lac->DrawLine(min_coin_c,ac1_adc,max_coin_c,ac1_adc);
 ccoin_ac->cd(2);
 hcoin_ac2->Draw("colz");
 lac->DrawLine(min_coin_c,ac2_adc,max_coin_c,ac2_adc); 


 TCanvas* cAC=new TCanvas("cAC","cAC");
 cAC->cd();
 ha1_a2->Draw

("colz");
 lac->DrawLine(min_ac1,ac2_adc,max_ac1,ac2_adc);
 lac->DrawLine(ac1_adc,min_ac2,ac1_adc,max_ac2);

 TCanvas* cmm=new TCanvas("cmm","cmm");
 cmm->cd();
 hmm->Draw();


 /*
   TCanvas* cAC1=new TCanvas("cAC1","cAC1");
   cAC1->Divide(1,3);
   cAC1->cd(1);
   grate_k_ac1->Draw("Ap");
   cAC1->cd(2);

   grate_pi_ac1->Draw("AP");
   cAC1->cd(3);
   grate_p_ac1->Draw("AP");
   m
   TCanvas* cAC2=new TCanvas("cAC2","cAC2");

   cAC2->Divide(1,3);
   cAC2->cd(1);
   grate_k_ac2->Draw("Ap");
   cAC2->cd(2);
   grate_pi_ac2->Draw("AP");
   cAC2->cd(3);
   grate_p_ac2->Draw("AP");
   

   TCanvas* cAC1_num=new TCanvas("cAC1_num","cAC1_num");
   cAC1_num->Divide(1,3);

   cAC1_num->cd(1);

   gsum_k_ac1->Draw("Ap");
   cAC1_num->cd(2);
   gsum_pi_ac1->Draw("AP");
   cAC1_num->cd(3);
   gsum_p_ac1->Draw("AP");
   
   TCanvas* cAC2_num=new TCanvas("cAC2_num","cAC2_num");
m
   cAC2_num->Divide(1,3);
   cAC2_num->cd(1);

   gsum_k_ac2->Draw("Ap");
   cAC2_num->cd(2);
   gsum_pi_ac2->Draw("AP");
   cAC2_num->cd(3);
   gsum_p_ac2->Draw("AP");
   

 */





 //================ Print Canvas =================================//
 
TString name;
 if(output_flag){

   // name.Form("./pdf/hdrogen_run_KID_eff_check.pdf");
   name.Form(ofname.c_str());
 cmm->Print(name+"[","pdf");
 cmm->Print(name,"pdf");
 cmm_cut->Print(name,"pdf");
 ccoin_ac->Print(name,"pdf");
 cAC->Print(name,"pdf");
 crpathl->Print(name,"pdf");
 crpathl->Print(name+"]","pdf");
 
}
 


 //====== Comment Out =====================//
 /*
 cout<<"======= Get Fit Parameters ============"<<endl;
 for(int i=0;i<iter_ac1;i++){
   cout<<"========="<<i<<"=================="<<endl;
 cout<<"AC1 threshold :"<<min_ac1+(ac1_adc-min_ac1)/iter_ac1*i<<endl;
 cout<<"Accidental BG p0: "<<p0_acc[i][0]<<endl;
 cout<<"Proton Fit parameters "<<endl;
 cout<<"mean "<<mean_p[i][0]<<endl;
 cout<<"sigma "<<sig_p[i][0]<<endl;


 cout<<"Pion Fit Parameters "<<endl;
 cout<<"mean "<<mean_pi[i][0]<<endl;
 cout<<"sigma "<<sig_pi[i][0]<<endl;
 cout<<"Kaon Fit Paramters "<<endl;
 cout<<"mean "<<mean_k[i][0]<<endl;
 cout<<"sigma "<<sig_k[i][0]<<endl;
 cout<<"Sum of Kaon "<<sum_k[i][0]<<endl;
 cout<<"Sum of Pion "<<sum_pi[i][0]<<endl;
 cout<<"Sum of Proton "<<sum_p[i][0]<<endl;
   cout<<"==========================="<<endl;
 }

 */



  cout<<"def_t3_p: "<<def_t1_p/tdc_time<<endl;
  cout<<Form(" : max p ac1[%d]: ",ac1_p)<<max_np[0]/tdc_time<<endl;
  cout<<Form(" : max p ac2[%d]: ",ac2_p)<<max_np[1]/tdc_time<<endl;
  cout<<"t1mean_p: "<<t1mean_p<<" : t2mean_p: "<<t2mean_p<<endl;

  cout<<"t1sig_p: "<<t1sig_p<<" : t2sig_p: "<<t2sig_p<<endl;
  cout<<"def_t3_pi: "<<def_t3_pi/tdc_time<<endl;
 cout<<Form(" : max pi ac1[%d]: ",ac1_pi)<<max_npi[0]/tdc_time<<endl;
 cout<<Form(" : max pi ac2[%d]: ",ac2_pi)<<max_npi[1]/tdc_time<<endl;
 cout<<"t1mean_pi: "<<t1mean_pi<<" : t2mean_pi: "<<t2mean_pi<<endl;
 cout<<"t1sig_pi: "<<t1sig_pi<<" : t2sig_pi: "<<t2sig_pi<<endl;
 
 cout<<"def_t3_k: "<<def_t3_k/tdc_time<<endl;
 cout<<Form(" : max k ac1[%d]: ",ac1_k)<<max_nk[0]/tdc_time<<endl;
 cout<<Form(" : max k ac2[%d]: ",ac2_k)<<max_nk[1]/tdc_time<<endl;

 cout<<"t1mean_k: "<<t1mean_k<<" : t2mean_k: "<<t2mean_k<<endl;
 cout<<"t1sig_k: "<<t1sig_k<<" : t2sig_k: "<<t2sig_k<<endl;
 cout<<"t3mean_k: "<<t3mean_k<<" :t3sig_k: "<<t3sig_k<<endl;


 gSystem->Exit(1);
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
