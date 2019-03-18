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
double s2f1_off(int i,char* ARM,int j,char* MODE,int KINE);
const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
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

  //------ Initial Parameters Setting ---------------//
  int ch; char* mode="H";
  int kine=1;// 1: hydrogen kinematics 2:tritium kinematics
  double tdc_time=58.0e-3;//[ns]
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  string ofroot = "rootfiles/hrs_test.root";

  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool root_flag=false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:r:bcop:GHT12"))!=-1){
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
  
    case 'r':
      ofroot = optarg;
      root_flag = true;
      break;
    
    case 'G':
    mode="G";
      break;
  
  case 'H':
    mode="H";
      break;

  case '1':
    tdc_time=56.23e-3;//[ns]
    kine=1;
      break;

  case '2':
    tdc_time=58e-3;//[ns]
    kine=2;
      break;
    
  case 'T':
    mode="T";
      break;
    
  
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
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


 if(draw_flag==0)gROOT->SetBatch(1);
 // gStyle->SetOptStat(000000000);
  //=============== ROOT File Mode ================//


 //TChain //
  TChain* T;
  if(mode=="G"){T=new TChain("tree"); }
  else {T=new TChain("T");}

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //     cout<<buf<<endl;
  }

 
 int evnt=T->GetEntries();
 cout<<"mode :"<<mode<<endl;
 cout<<"tdc_time : "<<tdc_time<<endl;
 cout<<"Get Entries: "<<evnt<<endl;



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
  double Rs2trpad[100],Ls2trpad[100];
  double ctime[100],DRT5;
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


 if(mode=="G"){
 T->SetBranchStatus("ctime",1);    
 T->SetBranchAddress("ctime",ctime);
 T->SetBranchStatus("DR.T5",1);    
 T->SetBranchAddress("DR.T5",&DRT5);

}

 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double min_coin,max_coin;
 double min_coin_c,max_coin_c;
 double ac1_adc,ac2_adc;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;

 // Range setting //

 if(mode=="H"){
   
   /*
    min_coin=-100.0; //test
    max_coin=100.0;//test
    min_coin_c=-50; //test
    max_coin_c=100.0;//test
   */
    min_coin=-10.0;
    max_coin=20.0;
    min_coin_c=-10;
    max_coin_c=20.0;
   

 min_ac1=0.0;
 max_ac1=5000.;
 min_ac2=0.0;
 max_ac2=20000.;
 min_adc=-500.0;
 max_adc=20000.;
//=== AC Threshold variable ===//

 ac1_adc=150.;
 ac2_adc=1000.;

 }else if(mode=="G"){
 min_coin=-20;
 max_coin=20.0;
 min_coin_c=-20;
 max_coin_c=20.0;
 min_ac1=0.0;
 max_ac1=30.;
 min_ac2=0.0;
 max_ac2=50.;
 min_adc=-5.0;
 max_adc=20.;
//==== AC Threshold variable (ACTH)===//
 ac1_adc=1.3;
 ac2_adc=13;
 
}


 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 
 int bin_beta=6000;
 int bin_adc=max_adc-min_adc;
 int bin_ac1=max_ac1-min_ac1; 
 int bin_ac2=max_ac2-min_ac2; 


 //=========== Coincidence Hist =================================//
 //-------- No correction -----------------//
 TH1F* hcoin=new TH1F("hcoin","Coincidence time S2R-S2L[ns] no  Parameter correction ",bin_coin,min_coin,max_coin);
 hcoin->SetTitle("w/o S2 No Correction Coincedence time R-S2 - L-S2; Coinc time [ns] ; Counts ");
 hcoin->SetTitleSize(0.1,"x");
 hcoin->SetTitleSize(0.1,"y");
 hcoin->SetLabelSize(0.05,"x");
 hcoin->SetLabelSize(0.05,"y");
 //-----------  Path Length correction --------------------//
 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] w/ Path length correction",bin_coin,min_coin,max_coin);
 hcoin_t->SetTitle("Coincedence time R-S2 - L-S2; Coinc time [ns] ; Counts ");
 hcoin_t->SetTitleSize(0.1,"x");
 hcoin_t->SetTitleSize(0.1,"y");
 hcoin_t->SetLabelSize(0.05,"x");
 hcoin_t->SetLabelSize(0.05,"y");
 //----------- Offset correction ------------------------//
 TH1F* hcoin_off=new TH1F("hcoin_off","Coincidence time S2R-S2L[ns] w/ S2 Offset Parameter correction "//,4000,-500.,0.);
			  ,bin_coin,min_coin,max_coin);
 hcoin_off->SetTitle("w/ S2 Offset Correction Coincedence time R-S2 - L-S2; Coinc time [ns] ; Counts ");
 hcoin_off->SetTitleSize(0.1,"x");
 hcoin_off->SetTitleSize(0.1,"y");
 hcoin_off->SetLabelSize(0.05,"x");
 hcoin_off->SetLabelSize(0.05,"y");
 //--------- after Correction -----------------------------//
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length && offset Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);

 //======================================================================//


 TH1F* hcoin_t1=new TH1F("hcoin_t1","Coincidence time S2R-S2L[ns] w/ AC1 cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t2=new TH1F("hcoin_t2","Coincidence time S2R-S2L[ns] w/ AC2 cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_t3=new TH1F("hcoin_t3","Coincidence time S2R-S2L[ns] w/ AC1 && AC2 cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_tk=new TH1F("hcoin_tk","Coincidence time S2R-S2L[ns] Kaoin cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_tp=new TH1F("hcoin_tp","Coincidence time S2R-S2L[ns] Proton cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_tpi=new TH1F("hcoin_pi","Coincidence time S2R-S2L[ns] Pion cut",bin_coin_c,min_coin_c,max_coin_c);
 double min_rpathl=28.5; double max_rpathl=29.5; int bin_rpathl=500;
 double min_lpathl=28.5; double max_lpathl=29.5; int bin_lpathl=500; 
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
 ha1_a2->SetTitle("AC1 vs AC2 ADC Hist;AC1 adc;AC2 adc; Counts");
 int bin_seg=16;
 double min_seg=0.0; double max_seg=15.;
 TH2F*hcoin_rseg=new TH2F("hcoin_rseg","Coinc Time vs R-S2seg Hist",bin_seg,min_seg,max_seg,bin_coin_c,min_coin_c,max_coin_c);

 TH2F*hcoin_lseg=new TH2F("hcoin_lseg","Coinc Time vs L-S2seg Hist",bin_seg,min_seg,max_seg,bin_coin_c,min_coin_c,max_coin_c);


 TH1F* hcoin_rs2[16]; 
 TH1F* hcoin_ls2[16];
 TH1F* hcoin_rs2_def[16]; 
 TH1F* hcoin_ls2_def[16];

 for(int i=0;i<16;i++){

   hcoin_rs2[i]=new TH1F(Form("hcoin_rs2[%d]",i),Form("Coincidence Time with R-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_ls2[i]=new TH1F(Form("hcoin_ls2[%d]",i),Form("Coincidence Time with L-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_rs2_def[i]=new TH1F(Form("hcoin_rs2_def[%d]",i),Form("Coincidence Time with R-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_ls2_def[i]=new TH1F(Form("hcoin_ls2_def[%d]",i),Form("Coincidence Time with L-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);


} 




 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 double mh;
 double m2; 
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rtof[16];
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off,Rs2_off2,Ls2_off2; 
 double Rs2_tcorr,Ls2_tcorr;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double tof_r,tof_l, tof_r2,tof_l2;
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr;
 double coin_offset;
 double coin,coin_tc2;
 double coin_off;
 double coin_path_offset;
 double coin_pc;

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

 //--- Coin Offset -----//
 double pathl_off,s2_offset;
 pathl_off=0.0;
 if(mode=="H"&&kine==2){
 coin_offset=498.;
 s2_offset=-489.0;
 pathl_off=-485.0; 

 

 }else if(mode=="H" && kine==1){
 pathl_off=-498.+30.-3.0+0.5;
 s2_offset=-500.0+25.;
 coin_offset=-41.35+498.;
}

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
 //rbeta=ppi/Epi; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 Rs2_off=s2f1_off(Rs2pads,"R",0,mode,kine);
 Ls2_off=s2f1_off(Ls2pads,"L",0,mode,kine);
 Rs2_off2=s2f1_off(Rs2pads,"R",1,mode,kine);
 Ls2_off2=s2f1_off(Ls2pads,"L",1,mode,kine);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 tof_r2=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off2)/2.0))*tdc_time;
 tof_l2=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off2)/2.0)*tdc_time;

 

 if(mode=="G"){
   coin_t=ctime[0];
   coin_tc=ctime[0];
 }else{
   /*
   tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
   tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
   tof_r2=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off2)/2.0))*tdc_time;
   tof_l2=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off2)/2.0)*tdc_time;
   */
   coin=(-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0*tdc_time-((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time-coin_offset; //coin no correct
   coin_pc=coin+rpath_corr-lpath_corr-pathl_off; // coin Path correct
   coin_t=tof_r-tof_l-coin_offset-s2_offset; // coin S2-Offset Correction
   coin_off=tof_r2-tof_l2-coin_offset; // coin s2 Offset Correct2
   coin_tc=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction
   coin_tc2=tof_r2-tof_l2+rpath_corr-lpath_corr-pathl_off-coin_offset; //coin Path & Offset2 Correction
   
  // coin_off=tof_r2-tof_l2+239.;//no offset version
  
}


//====== Cut condition ========================// 
   cut_ac1=false;
   cut_ac2=false;
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   coin_trig=false;
   cut_track=false;
   cut_s0=false;
   // if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
   if(Ra1sum>ac1_adc)cut_ac1=true;
   if(Ra2sum<ac2_adc)cut_ac2=true;
   //  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
   // if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
   if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true;
   // if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
   // if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
   //  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;


  if(mode=="G"){
     if(ctime[0]>-100)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true;
   }else{

 if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
 if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
 if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
 if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
 if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
 if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
 if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true;

   
}



 //=======================================//



   //==========================================//
   //========= Fill Hist =====================//
   //========================================//



 if(coin_trig && cut_vz && cut_track){
 hcoin_rpathl->Fill(rpathl,coin_t);
 hcoin_rpathl_cc->Fill(rpathl,coin_tc);
 hcoin_lpathl->Fill(lpathl,coin_t);
 hcoin_lpathl_cc->Fill(lpathl,coin_tc);
 }
  
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){
 hcoin->Fill(coin);
 hcoin_t->Fill(coin_t);
 hcoin_off->Fill(coin_off);
 hcoin_tc->Fill(coin_tc);
 hcoin_ac1->Fill(coin_tc,Ra1sum);
 hcoin_ac2->Fill(coin_tc,Ra2sum);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2

 hcoin_rseg->Fill(Rs2pads,coin_tc);
 hcoin_lseg->Fill(Ls2pads,coin_tc);

 if(Ra1sum>ac1_adc && Ra2sum>ac2_adc){
 hcoin_tpi->Fill(coin_tc); 
 hcoin_rs2[Rs2pads]->Fill(coin_tc); 
 hcoin_ls2[Ls2pads]->Fill(coin_tc);
 // hcoin_rs2_def[Rs2pads]->Fill(coin_off); //defolt
 // hcoin_ls2_def[Ls2pads]->Fill(coin_off); //defolt
 hcoin_rs2_def[Rs2pads]->Fill(coin); //defolt
 hcoin_ls2_def[Ls2pads]->Fill(coin); //defolt


 }
 if(Ra1sum<ac1_adc &&Ra2sum>ac2_adc)hcoin_tk->Fill(coin_tc); // Kaon Cut

 if(Ra1sum<ac1_adc && Ra2sum<ac2_adc)hcoin_tp->Fill(coin_tc); 
 }

 if(coin_trig && cut_ac1 && cut_vz && cut_lpathl && cut_rpathl && cut_track){
  hcoin_t2->Fill(coin_tc); //ac1 cut fix
 
 }
 if(coin_trig && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track){
  hcoin_t1->Fill(coin_tc); // ac2 cut fix
 

 }
 if(coin_trig  && cut_ac1 && cut_ac2  && cut_vz && cut_lpathl && cut_rpathl && cut_track){
   hcoin_t3->Fill(coin_tc); //ac1 & ac2 cut

        }

 

 }





 //======= Kaon Event Analysis ============================//                   


 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_sig_k,def_mean_k,def_acc;


 if(mode=="H" && kine==2){
 def_sig_p=0.852; def_mean_p=-13.82;
 def_sig_pi=0.443; def_mean_pi=-2.63;
 def_sig_k=0.644; def_mean_k=-5.55;
 def_acc=27.7;
 }else if(mode=="H" && kine==1){
def_sig_p=0.852; def_mean_p=-13.82;
 def_sig_pi=0.443; def_mean_pi=-2.63;
 def_sig_k=0.644; def_mean_k=-5.55;
 def_acc=27.7;
 }else if(mode=="G"){
 def_sig_p=0.852; def_mean_p=11.06;
 def_sig_pi=0.4; def_mean_pi=0.0;
 def_sig_k=0.644; def_mean_k=3.16;
 def_acc=22.7;
}






 TH1F* hcoin_tx=(TH1F*)hcoin_tk->Clone();
 hcoin_tx->SetName("hcoin_tx");

 TF1* facc=new TF1("facc","[0]",min_coin_c,max_coin_c);
 facc->SetNpx(2000);
 hcoin_tx->Fit("facc","R","",min_coin_c,min_coin_c+3.0);
 double p0_acc=facc->GetParameter(0);
 // double p1_acc=facc->GetParameter(1);                                        
 TF1* fp=new TF1("fp","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp->SetNpx(2000);
 fp->FixParameter(3,p0_acc);
 fp->SetParameter(1,0.0);
 fp->SetParameter(2,0.788);
 hcoin_tx->Fit("fp","R","",-5.,5);
 double n_p=fp->GetParameter(0);
 double mean_p=fp->GetParameter(1);
 double sig_p=fp->GetParameter(2);
 TF1* fpi=new TF1("fpi","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi->SetNpx(2000);
 fpi->SetParameter(3,p0_acc);
 fpi->SetParameter(1,11.);
 fpi->SetParameter(2,0.384);
 hcoin_tx->Fit("fpi","R","",10.,15);
 double n_pi=fpi->GetParameter(0);
 double mean_pi=fpi->GetParameter(1);
 double sig_pi=fpi->GetParameter(2);
 TF1* fk=new TF1("fk","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk->SetNpx(2000);
 fk->SetParameter(3,p0_acc);
 fk->SetParameter(1,9.5);
 //fk->SetParLimits(1,5.4,5.5);                                                 
 fk->SetParameter(2,0.844);
 hcoin_tx->Fit("fk","Rb","",5.0,10.0);
 double n_k=fk->GetParameter(0);
 double mean_k=fk->GetParameter(1);
 double sig_k=fk->GetParameter(2);


 TF1* fcoin=new TF1("fcoin","gausn(0)+gausn(3)+gausn(6)+pol0(9)");
 fcoin->SetNpx(2000);
 fcoin->SetParameters(n_pi,mean_pi,sig_pi,n_k,mean_k,sig_k,n_p,mean_p,sig_p,p0_acc);
 fcoin->SetParLimits(4,mean_k-0.5*sig_k,mean_k+0.5*sig_k);
 fcoin->SetParLimits(5,0.8*sig_k,1.2*sig_k);
 hcoin_tx->Fit("fcoin","","",min_coin_c,max_coin_c);

 n_pi=fcoin->GetParameter(0); mean_pi=fcoin->GetParameter(1);  
 sig_pi=fcoin->GetParameter(2);
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

 n=hcoin_t->GetEntries();
 nac1=hcoin_t1->GetEntries();
 nac2=hcoin_t2->GetEntries();
 nac3=hcoin_t3->GetEntries();






 //============ COINCIDENCE SEGMENT CORRECTION ================//

 double min_pi,max_pi,rn_pi[16],rsig_pi[16],rmean_pi[16],ln_pi[16],lsig_pi[16],lmean_pi[16],rmean_pi_err[16],lmean_pi_err[16];
 double diff_r[16],diff_l[16],r_pmax[16],l_pmax[16],nl_pmax[16],nr_pmax[16],bin_rmax[16],bin_lmax[16],bin_rxmax[16],bin_rymax[16],bin_lxmax[16],bin_lymax[16];
 double rn_pi_def[16],rsig_pi_def[16],rmean_pi_def[16],ln_pi_def[16],lsig_pi_def[16],lmean_pi_def[16],rmean_pi_err_def[16],lmean_pi_err_def[16];
 double diff_r_def[16],diff_l_def[16],r_pmax_def[16],l_pmax_def[16],nl_pmax_def[16],nr_pmax_def[16],bin_rmax_def[16],bin_lmax_def[16],bin_rxmax_def[16],bin_rymax_def[16],bin_lxmax_def[16],bin_lymax_def[16];
 min_pi=0., max_pi=25;
 TF1* ffit_r[16];
 TF1* ffit_l[16];
 TF1* ffit_r_def[16];
 TF1* ffit_l_def[16];


 TGraphErrors* gs2r_mean=new TGraphErrors();
 gs2r_mean->SetName("gs2r_mean");
 gs2r_mean->SetTitle("coincidence tiem of r-s2 segmments; segment num ;coin time [ns]");
 gs2r_mean->GetXaxis()->SetTitleSize(0.05);
 gs2r_mean->GetYaxis()->SetTitleSize(0.05);
 gs2r_mean->GetXaxis()->SetLabelSize(0.05);
 gs2r_mean->GetYaxis()->SetLabelSize(0.04);
 TGraphErrors* gs2l_mean=new TGraphErrors();
 gs2l_mean->SetName("gs2l_mean");
 gs2l_mean->SetTitle("coincidence tiem of l-s2 segments dependence ; segment num ;coin time [ns] ");
 gs2l_mean->GetXaxis()->SetTitleSize(0.05);
 gs2l_mean->GetYaxis()->SetTitleSize(0.05);
 gs2l_mean->GetXaxis()->SetLabelSize(0.05);
 gs2l_mean->GetYaxis()->SetLabelSize(0.04);
TGraphErrors* gs2r_mean_def=new TGraphErrors();
TGraphErrors* gs2l_mean_def=new TGraphErrors();

 gs2r_mean->GetYaxis()->SetRangeUser(-2.0,2.0);
 gs2l_mean->GetYaxis()->SetRangeUser(-2.0,2.0);


 //=======================================================//



ffit_r[7] =new TF1(Form("ffit_r[%d]",7),"gaus",min_coin_c,max_coin_c);
ffit_l[7] =new TF1(Form("ffit_l[%d]",7),"gaus",min_coin_c,max_coin_c);
ffit_r_def[7] =new TF1(Form("ffit_r_def[%d]",7),"gaus",min_coin_c,max_coin_c);
ffit_l_def[7] =new TF1(Form("ffit_l_def[%d]",7),"gaus",min_coin_c,max_coin_c);

 hcoin_rs2[7]->Fit(Form("ffit_r[%d]",7),"0","",min_pi,max_pi);
 rmean_pi[7]=ffit_r[7]->GetParameter(1);
 rmean_pi_err[7]=ffit_r[7]->GetParError(1);
 hcoin_ls2[7]->Fit(Form("ffit_l[%d]",7),"0","",min_pi,max_pi);
 lmean_pi[7]=ffit_l[7]->GetParameter(1);
 lmean_pi_err[7]=ffit_l[7]->GetParError(1);


 hcoin_rs2_def[7]->Fit(Form("ffit_r_def[%d]",7),"0","",min_pi,max_pi);
 rmean_pi_def[7]=ffit_r_def[7]->GetParameter(1);
 rmean_pi_err_def[7]=ffit_r_def[7]->GetParError(1);
 hcoin_ls2_def[7]->Fit(Form("ffit_l_def[%d]",7),"0","",min_pi,max_pi);
 lmean_pi_def[7]=ffit_l_def[7]->GetParameter(1);
 lmean_pi_err_def[7]=ffit_l_def[7]->GetParError(1);

 //=======================================================//


   double rs2_off[16],ls2_off[16];
   double rs2_off_def[16],ls2_off_def[16];
   for(int i=1;i<15;i++){
   
   // for(int j=0;j<2;j++){
   ffit_r[i] =new TF1(Form("ffit_r[%d]",i),"gaus",min_coin_c,max_coin_c);
   ffit_l[i] =new TF1(Form("ffit_l[%d]",i),"gaus",min_coin_c,max_coin_c);
   ffit_r_def[i] =new TF1(Form("ffit_r_def[%d]",i),"gaus",min_coin_c,max_coin_c);
   ffit_l_def[i] =new TF1(Form("ffit_l_def[%d]",i),"gaus",min_coin_c,max_coin_c);


   //----------- l-hrs -------------------------//

   bin_rmax[i]= hcoin_rs2[i]->GetMaximumBin();
   r_pmax[i]= hcoin_rs2[i]->GetXaxis()->GetBinCenter(bin_rmax[i]);
   nr_pmax[i]= hcoin_rs2[i]->GetBinContent(bin_rmax[i]);
   rsig_pi[i]=3.5;   

   bin_rmax_def[i]= hcoin_rs2_def[i]->GetMaximumBin();
   r_pmax_def[i]= hcoin_rs2_def[i]->GetXaxis()->GetBinCenter(bin_rmax_def[i]);
   nr_pmax_def[i]= hcoin_rs2_def[i]->GetBinContent(bin_rmax_def[i]);
   rsig_pi_def[i]=3.5;   

   ffit_r[i]->SetParameter(0,nr_pmax[i]);
   ffit_r[i]->SetParameter(1,r_pmax[i]);
   ffit_r[i]->FixParameter(2,rsig_pi[i]);
   ffit_r[i]->SetParLimits(0,0.8*nr_pmax[i],1.2*nr_pmax[i]);
   ffit_r[i]->SetParLimits(1,r_pmax[i]-1.0,r_pmax[i]+1.0);
   ffit_r[i]->SetParLimits(2,0.5,4.0);

   ffit_r_def[i]->SetParameter(0,nr_pmax_def[i]);
   ffit_r_def[i]->SetParameter(1,r_pmax_def[i]);
   ffit_r_def[i]->FixParameter(2,rsig_pi_def[i]);
   ffit_r_def[i]->SetParLimits(0,0.8*nr_pmax_def[i],1.2*nr_pmax_def[i]);
   ffit_r_def[i]->SetParLimits(1,r_pmax_def[i]-1.0,r_pmax_def[i]+1.0);
   ffit_r_def[i]->SetParLimits(2,0.5,4.0);
    
   hcoin_rs2[i]->Fit(Form("ffit_r[%d]",i),"Rb","",min_pi,max_pi);
   rmean_pi[i]=ffit_r[i]->GetParameter(1);
   rmean_pi_err[i]=ffit_r[i]->GetParError(1);
   nr_pmax[i]=ffit_r[i]->GetParameter(0);
   rsig_pi[i]=ffit_r[i]->GetParameter(2);
   nl_pmax[i]=hcoin_ls2[i]->GetYaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());


   hcoin_rs2_def[i]->Fit(Form("ffit_r_def[%d]",i),"Rb","",min_pi,max_pi);
   rmean_pi_def[i]=ffit_r_def[i]->GetParameter(1);
   rmean_pi_err_def[i]=ffit_r_def[i]->GetParError(1);
   nr_pmax_def[i]=ffit_r_def[i]->GetParameter(0);
   rsig_pi_def[i]=ffit_r_def[i]->GetParameter(2);

   //----------- l-hrs -------------------------//

   bin_lmax[i]= hcoin_ls2[i]->GetMaximumBin();
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(bin_lmax[i]);
   nl_pmax[i]= hcoin_ls2[i]->GetBinContent(bin_lmax[i]);
   lsig_pi[i]=3.5;   

   bin_lmax_def[i]= hcoin_ls2_def[i]->GetMaximumBin();
   l_pmax_def[i]= hcoin_ls2_def[i]->GetXaxis()->GetBinCenter(bin_lmax_def[i]);
   nl_pmax_def[i]= hcoin_ls2_def[i]->GetBinContent(bin_lmax_def[i]);
   lsig_pi_def[i]=3.5;   



   ffit_l[i]->SetParameter(0,nl_pmax[i]);
   ffit_l[i]->SetParameter(1,l_pmax[i]);
   ffit_l[i]->FixParameter(2,lsig_pi[i]);
   ffit_l[i]->SetParLimits(0,0.8*nl_pmax[i],1.2*nl_pmax[i]);
   ffit_l[i]->SetParLimits(1,l_pmax[i]-1.0,l_pmax[i]+1.0);
   ffit_l[i]->SetParLimits(2,0.5,4.0);

   ffit_l_def[i]->SetParameter(0,nl_pmax_def[i]);
   ffit_l_def[i]->SetParameter(1,l_pmax_def[i]);
   ffit_l_def[i]->FixParameter(2,lsig_pi_def[i]);
   ffit_l_def[i]->SetParLimits(0,0.8*nl_pmax_def[i],1.2*nl_pmax_def[i]);
   ffit_l_def[i]->SetParLimits(1,r_pmax_def[i]-1.0,r_pmax_def[i]+1.0);
   ffit_l_def[i]->SetParLimits(2,0.5,4.0);
   

   hcoin_ls2[i]->Fit(Form("ffit_l[%d]",i),"Rb","",min_pi,max_pi);
   lmean_pi[i]=ffit_l[i]->GetParameter(1);
   lmean_pi_err[i]=ffit_l[i]->GetParError(1);
   nl_pmax[i]=ffit_l[i]->GetParameter(0);
   lsig_pi[i]=ffit_l[i]->GetParameter(2);
   nl_pmax[i]=hcoin_ls2[i]->GetYaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());


   hcoin_ls2_def[i]->Fit(Form("ffit_l_def[%d]",i),"Rb","",min_pi,max_pi);
   lmean_pi_def[i]=ffit_l_def[i]->GetParameter(1);
   lmean_pi_err_def[i]=ffit_l_def[i]->GetParError(1);
   nl_pmax_def[i]=ffit_l_def[i]->GetParameter(0);
   lsig_pi_def[i]=ffit_l_def[i]->GetParameter(2);
   // lmean_pi_err[i]=ffit_l[i]->Getparerror(1);


 gs2r_mean->SetPoint(i,i,(rmean_pi[i]-rmean_pi[7]));
 gs2r_mean_def->SetPoint(i,i,(rmean_pi_def[i]-rmean_pi_def[7]));
// gs2r_mean->setpointerror(i,0,rmean_pi_err[i]);

 gs2l_mean->SetPoint(i,i,lmean_pi[i]-lmean_pi[7]);
 gs2l_mean_def->SetPoint(i,i,lmean_pi_def[i]-lmean_pi_def[7]);
 // gs2l_mean->setpointerror(i,0,lmean_pi_err[i]);

 diff_r[i]=-(rmean_pi[i]-rmean_pi[7])/tdc_time*2;
 diff_l[i]=-(lmean_pi[i]-lmean_pi[7])/tdc_time*2;
 rs2_off[i]=-(rmean_pi[i]-rmean_pi[7])/tdc_time*2+s2f1_off(i,"R",1,mode,kine);
 ls2_off[i]=(lmean_pi[i]-lmean_pi[7])/tdc_time*2+s2f1_off(i,"L",1,mode,kine);
 diff_r_def[i]=-(rmean_pi_def[i]-rmean_pi_def[7])/tdc_time*2;
 diff_l_def[i]=-(lmean_pi_def[i]-lmean_pi_def[7])/tdc_time*2;
 rs2_off_def[i]=-(rmean_pi_def[i]-rmean_pi_def[7])/tdc_time*2+s2f1_off(i,"R",0,mode,kine);
 ls2_off_def[i]=(lmean_pi_def[i]-lmean_pi_def[7])/tdc_time*2+s2f1_off(i,"L",0,mode,kine);


}
 // }



 //=============================================================//





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
 // hcoin_ac1->Draw();
 lac->DrawLine(min_coin_c,ac1_adc,max_coin_c,ac1_adc);
 ccoin_ac->cd(2);
 hcoin_ac2->Draw("colz");
 // hcoin_ac2->Draw();
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




 TCanvas* cs2seg=new TCanvas("cs2seg","s2 seg vs coincidence time graph rhrs[1] && lhrs[2] blue is w/o offSet black is w/ offSet ");
 cs2seg->Divide(2,1);
 cs2seg->cd(1);
 gs2r_mean->SetMarkerStyle(21);
 gs2r_mean->SetMarkerColor(kBlue);
 gs2r_mean->SetMarkerSize(1.0);
 gs2r_mean->Draw("ap");
 cs2seg->cd(2);
 gs2l_mean->SetMarkerStyle(21);
 gs2l_mean->SetMarkerColor(kBlue);
 gs2l_mean->SetMarkerSize(1.0); 
 gs2l_mean->Draw("ap");
 cs2seg->cd(1);
 gs2r_mean_def->SetMarkerStyle(21);
 gs2r_mean_def->SetMarkerColor(1);
 gs2r_mean_def->SetMarkerSize(1.0);
 gs2r_mean_def->Draw("p");
 cs2seg->cd(2);
 gs2l_mean_def->SetMarkerStyle(21);
 gs2l_mean_def->SetMarkerColor(1);
 gs2l_mean_def->SetMarkerSize(1.0); 
 gs2l_mean_def->Draw("p");

 

 TCanvas* crs2_pi_def=new TCanvas("crs2_pi_def","crs2_pi_def");
 crs2_pi_def->Divide(4,4);
 TCanvas* cls2_pi_def=new TCanvas("cls2_pi_def","cls2_pi_def");
 cls2_pi_def->Divide(4,4);
 TCanvas* crs2_pi=new TCanvas("crs2_pi","crs2_pi");
 crs2_pi->Divide(4,4);
 TCanvas* cls2_pi=new TCanvas("cls2_pi","cls2_pi");
 cls2_pi->Divide(4,4);
 for(int i=0;i<16;i++){
 crs2_pi->cd(i+1);
 hcoin_rs2[i]->SetLineColor(kBlue);
 hcoin_rs2[i]->Draw();
 crs2_pi_def->cd(i+1);
 hcoin_rs2_def[i]->SetLineColor(1);
 hcoin_rs2_def[i]->Draw();
 cls2_pi->cd(i+1);
 hcoin_ls2[i]->SetLineColor(kBlue);
 hcoin_ls2[i]->Draw();
 cls2_pi_def->cd(i+1);
 hcoin_ls2_def[i]->SetLineColor(1);
 hcoin_ls2_def[i]->Draw();
}



 TCanvas* cseg=new TCanvas("cseg","cseg");
 cseg->Divide(1,2);
 cseg->cd(1); 
 hcoin_rseg->Draw("colz");
 cseg->cd(2); 
 hcoin_lseg->Draw("colz");


 
 /*
 TCanvas* cseg2=new TCanvas("cseg2","cseg2");
 cseg2->Divide(1,2);
 cseg2->cd(1);
 hcoin_rs2[8]->SetLineColor(9);
 hcoin_rs2[8]->Draw();
 cseg2->cd(2);
 hcoin_ls2[8]->SetLineColor(9);
 hcoin_ls2[8]->Draw();
 for(int i=0;i<16;i++){
   if(i==8){
   }else{
 cseg2->cd(1);
 hcoin_rs2[i]->SetLineColor(i+1);
 hcoin_rs2[i]->Draw("same");
 cseg2->cd(2);
 hcoin_ls2[i]->SetLineColor(i+1);
 hcoin_ls2[i]->Draw("same");
}
 }
 
 */
 
 /*
 TCanvas *cnum=new TCanvas("cnum","cnum");
 cnum->cd();
 fcoin->SetLineColor(kred);
 fk->SetLineColor(4);
 fk->SetFillColor(4);
 fk->SetFillStyle(3001);
 fpi->SetFillStyle(3001);
 fpi->SetLineColor(46);
 fpi->SetFillColor(46);
 fp->SetFillStyle(3001);
 fp->SetLineColor(8);
 fp->SetFillColor(8);
 hcoin_tx->Draw();
 fcoin->Draw("same");
 fk->Draw("same");
 fp->Draw("same");
 fpi->Draw("same");
 
*/
 
 TCanvas* cac=new TCanvas("cac","cac");
 cac->cd();
 ha1_a2->Draw("colz");
 //  ha1_a2->Draw();
 lac->DrawLine(min_ac1,ac2_adc,max_ac1,ac2_adc);
 lac->DrawLine(ac1_adc,min_ac2,ac1_adc,max_ac2);

 TCanvas* c_pathc =new TCanvas("c_pathc","c_pathc");
 c_pathc->Divide(1,2);
 c_pathc->cd(1);
 hcoin_tc->Draw();
 c_pathc->cd(2);
 hcoin->Draw();

 TCanvas *crpathl=new TCanvas("crpathl","crpathl");
 crpathl->Divide(2,2);
 crpathl->cd(1);
 hcoin_rpathl->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin,rpathl_cutmin,max_coin); 
 lac->DrawLine(rpathl_cutmax,min_coin,rpathl_cutmax,max_coin);
 crpathl->cd(2);
 hcoin_lpathl->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin,lpathl_cutmin,max_coin); 
 lac->DrawLine(lpathl_cutmax,min_coin,lpathl_cutmax,max_coin); 
 crpathl->cd(3);
 hcoin_rpathl_cc->Draw("colz");
 lac->DrawLine(rpathl_cutmin,min_coin_c,rpathl_cutmin,max_coin_c); 
 lac->DrawLine(rpathl_cutmax,min_coin_c,rpathl_cutmax,max_coin_c);
 crpathl->cd(4);
 hcoin_lpathl_cc ->Draw("colz");
 lac->DrawLine(lpathl_cutmin,min_coin_c,lpathl_cutmin,max_coin_c); 
 lac->DrawLine(lpathl_cutmax,min_coin_c,lpathl_cutmax,max_coin_c); 

 TCanvas * ck=new TCanvas("ck");
 ck->Divide(2,1);
 ck->cd(1);
 hcoin_t2->Draw();
 ck->cd(2);
 hcoin_t1->Draw();

 TCanvas* ccoin0=new TCanvas("ccoin0","coincidence time hist w/o & w/ offSet correction");

 ccoin0->cd();
 hcoin_tc->SetLineColor(4);
 hcoin_tc->SetFillColor(4);
 hcoin_tc->SetFillStyle(3001);
 hcoin_tc->Draw(); 

 hcoin->SetFillStyle(3001);
 hcoin->SetLineColor(8);
 hcoin->SetFillColor(8);
 hcoin->Draw("same");


 TCanvas* ccoin_corr=new TCanvas("ccoin_corr","coincidence time correction hist");
 ccoin_corr->Divide(2,2);
 //hcoin_tc->Draw();
 //---- no correction -----// 
 ccoin_corr->cd(1);
 hcoin->SetFillStyle(3001);
 hcoin->SetLineColor(8);
 hcoin->SetFillColor(8);
 hcoin->Draw();
 //---- s2 offSet correction ---//
 ccoin_corr->cd(2);
 hcoin_off->SetFillStyle(3001);
 hcoin_off->SetLineColor(2);
 hcoin_off->SetFillColor(2);
 hcoin_off->Draw();

//---- path length correction ---//
 ccoin_corr->cd(3);
 hcoin_t->SetFillStyle(3001);
 hcoin_t->SetLineColor(3);
 hcoin_t->SetFillColor(3);
 hcoin_t->Draw();

 //----- path & offSet correction ----//
 ccoin_corr->cd(4);
 hcoin_tc->SetLineColor(4);
 hcoin_tc->SetFillColor(4);
 hcoin_tc->SetFillStyle(3001);
 hcoin_tc->Draw();
 

TCanvas* c0=new TCanvas("c0","c0");

 c0->cd();
 fk->SetLineColor(4);
 fk->SetFillColor(4);
 fk->SetFillStyle(3001);
 fpi->SetFillStyle(3001);
 fpi->SetLineColor(46);
 fpi->SetFillColor(46);
 fp->SetFillStyle(3001);
 fp->SetLineColor(8);
 fp->SetFillColor(8);
 hcoin_tx->Draw();
 fcoin->Draw("same");
 fk->Draw("same");
 fp->Draw("same");
 fpi->Draw("same");


 //================ print canvas =================================//
 
 
 TString name;
 if(output_flag){
 name.Form(ofname.c_str());
 ccoin_ac->Print(name+"[","pdf");
 ccoin_ac->Print(name,"pdf");
 ccoin0->Print(name,"pdf");
 ccoin_corr->Print(name,"pdf");
 cac->Print(name,"pdf");
 cseg->Print(name,"pdf");
 // cseg2->Print(name,"pdf");
 c0->Print(name,"pdf");
 // cnum->Print(name,"pdf");
 cs2seg->Print(name,"pdf");
 crs2_pi_def->Print(name,"pdf");
 crs2_pi->Print(name,"pdf"); 
 cls2_pi_def->Print(name,"pdf");
 cls2_pi->Print(name,"pdf");
 c_pathc->Print(name,"pdf");
 crpathl->Print(name,"pdf");
 crpathl->Print(name +"]","pdf");
}




 //====== comment out =====================//

 cout<<" event w/o ac cut : "<<n<<endl;
 cout<<" event w/ ac1 cut : "<<nac1<<endl;
 cout<<" event w/ ac2 cut : "<<nac2<<endl;
 cout<<" event w/ ac1 & ac2 cut : "<<nac3<<endl;
 
 cout<<"======= get fit parameters ============"<<endl;
 cout<<"accidental bg p0: "<<p0_acc<<endl;
 cout<<"proton fit parameters "<<endl;
 cout<<"mean "<<mean_p<<endl;
 cout<<"sigma "<<sig_p<<endl;
 cout<<"pion fit parameters "<<endl;
 cout<<"mean "<<mean_pi<<endl;
 cout<<"sigma "<<sig_pi<<endl;
 cout<<"kaon fit paramters "<<endl;
 cout<<"mean "<<mean_k<<endl;
 cout<<"sigma "<<sig_k<<endl;
 cout<<"sum of kaon "<<sum_k<<endl;
 //cout<<"sum of kaon in hcoin_tk "<<sum_k2<<endl;
 cout<<"sum of pion "<<sum_pi<<endl;
 cout<<"sum of proton "<<sum_p<<endl;
 
 for(int i=0;i<16;i++){
 cout<<Form("rmin_pi[%d]",i)<<rmean_pi[i]<<endl;
 cout<<Form("lmin_pi[%d]",i)<<lmean_pi[i]<<endl;
 cout<<Form("rmin_pi_diff[%d]-[7]",i)<<rmean_pi[i]-rmean_pi[7]<<endl;
 cout<<Form("lmin_pi_diff[%d]-[7]",i)<<lmean_pi[i]-lmean_pi[7]<<endl;
 cout<<Form("rs2[%d] t0 diff [ch]: ",i)<<diff_r[i]<<endl;
 cout<<Form("ls2[%d] t0 diff [ch]: ",i)<<diff_l[i]<<endl;
 cout<<Form("rs2[%d] [ch]: ",i)<<rs2_off[i]<<endl;
 cout<<Form("ls2[%d] [ch]: ",i)<<ls2_off[i]<<endl;
 cout<<Form("r_pmax[%d]",i)<<r_pmax[i]<<endl; 
 cout<<Form("l_pmax[%d]",i)<<l_pmax[i]<<endl;
 cout<<Form("nr_pmax[%d]",i)<<nr_pmax[i]<<endl;
 cout<<Form("rsig_pi[%d]",i)<< rsig_pi[i]<<endl;

 }

 cout<<"====== parameters set ================="<<endl;
  cout<<"===== Rs2_off ======"<<endl;
 for(int i=0;i<16;i++){
  
   cout<<rs2_off[i]<<",";
}
 cout<<endl;
  cout<<"===== Ls2_off ======"<<endl;
for(int i=0;i<16;i++){
  
   cout<<ls2_off[i]<<",";
}
 cout<<endl; 


 if(root_flag==1){
 TFile * ofp=new TFile(Form("%s",ofroot.c_str()),"recreate");
 TTree* T=new TTree("T",Form("Hydrogen run root Phase %d",kine));
 T->Draw("counts>>hcoin_ac1");
 T->Draw("counts>>hcoin_ac2");
 T->Draw("counts>>ha1_a2");
 T->Draw("counts>>hcoin");
 T->Draw("counts>>hcoin_t");
 T->Draw("counts>>hcoin_t2");
 T->Draw("counts>>hcoin_t3");
 T->Draw("counts>>hcoin_tc");
 T->Draw("counts>>hcoin_tx");
 T->Write();
 ofp->Close();
}



 if(draw_flag==0)gSystem->Exit(1);
 theApp->Run();
 return 0;



}


//==============================================//
//========== defined function ==================//
//=============================================//
double s2f1_off(int i,char* ARM,int j,char* MODE,int KINE){



  double Ls2_off[16]; 
  double Rs2_off[16];
 


  if(j==0){
    if(MODE=="H" && KINE==1){

    
 double  RS2_off_H1[16]={-16911.4,-16864.9,-16900,-16897.6,-16874.8,-16869.3,-16901.1,-16876.8,-16895.6,-16860.3,-16892.6,-16885,-16847.3,-16843.3,-16838.4,-16882.6};
 double  LS2_off_H1[16]={-25336.9,-25385.7,-25367,-25392.2,-25391,-25386.3,-25422,-25428.9,-25415.2,-25425,-25438,-25381,-25394.4,-25417.5,-25432.8,-26082.1};
    

      //---------------- Pion Offset ------------------//
      /*
 double  RS2_off_H1[16]={-16911.4,-16853.7,-16893.9,-16893.7,-16871.5,-16868,-16900.5,-16876.8,-16897.9,-16865.4,-16898.9,-16892.5,-16855.7,-16856.5,-16851.7,-16882.6};
 double  LS2_off_H1[16]={-25336.9,-25391.4,-25375.8,-25397.7,-25391.7,-25387.3,-25422,-25428.9,-25423,-25431.9,-25445,-25386.3,-25381.2,-25413,-25426.5,-26082.1};
      */
     //-----------------------------------------------//
  Ls2_off[i]=LS2_off_H1[i];
  Rs2_off[i]=RS2_off_H1[i];
    }
 else  if(MODE=="H" && KINE==2){
   /* 
 double  RS2_off_H2[16]={-16911.4,-16864.9,-16900,-16897.6,-16874.8,-16869.3,-16901.1,-16876.8,-16895.6,-16860.3,-16892.6,-16885,-16847.3,-16843.3,-16838.4,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25385.7,-25367,-25392.2,-25391,-25386.3,-25422,-25428.9,-25415.2,-25425,-25438,-25381,-25394.4,-25417.5,-25432.8,-26082.1};
   */
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  Ls2_off[i]=LS2_off_H2[i];
  Rs2_off[i]=RS2_off_H2[i];
    }

  }else if(j==1){

 double  RS2_off_a[16]={-16911.4,-16861.7,-16900.5,-16899,-16875.6,-16870.3,-16901.2,-16876.8,-16896.5,-16862,-16894.3,-16886.1,-16847.9,-16846.9,-16840.3,-16882.6};
 double  LS2_off_a[16]={-25336.9,-25391,-25375.9,-25397.9,-25392.1,-25387.9,-25422.3,-25428.9,-25422.5,-25431.8,-25444.2,-25385.8,-25381.4,-25409.1,-25426.6,-26082.1};
  
  Ls2_off[i]=LS2_off_a[i];
  Rs2_off[i]=RS2_off_a[i];

}
  


 

 

  double s2f1_offset;
  if(ARM=="R")s2f1_offset=Rs2_off[i];
  else  if(ARM=="L")s2f1_offset=Ls2_off[i];
  else {cout<<"false read out !!"<<endl;}
  return s2f1_offset;
}




