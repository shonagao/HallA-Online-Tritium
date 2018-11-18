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
const  int nth=3; //th num
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




//===================================================================//
//============================= Main ================================//
//===================================================================//

int main(int argc, char** argv){

  int ch; char* mode;
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:GHT"))!=-1){
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
  
    case 'G':
    mode="G";
      break;
  
    case 'H':
    mode="H";
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



  //=============== ROOT File Mode ================//
  /*if(ifname.c_str()=="../run_list/coin_H2_1.root")mode="G";
  else if(ifname.c_str()=="../run_list/Lambda_small.list" ||ifname.c_str()=="../run_list/Lambda_test.list")mode="H_1";
  else{cout<<"false to read mode types ";};
  */


 //double tdc_time=56.23e-3;//[ns]
 double tdc_time=58.e-3;//[ns] 

/*
 if(mode=="G"){tdc_time=56.23e-3;}
 else if(mode=="H"){tdc_time=56.23e-3;}
 else if(mode=="T"){tdc_time=58e-3;}
 */


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
    istringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //  cout<<buf<<endl;
  }


  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;
  int evnt=T->GetEntries();
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
  //---- Gogami root ---------//
  double ctime[1000];
  double DRT5;


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

 double ac1_adc[nth],ac2_adc[nth];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;

 if(mode=="H" || mode=="T"){
 min_coin=-20;
 max_coin=0.0;
 // min_coin_c=-20;
 // max_coin_c=0.0;
  min_coin_c=-100;
  max_coin_c=100.0;
 min_ac1=0.0;
 max_ac1=5000.;
 min_ac2=0.0;
 max_ac2=20000.;
 min_adc=-500.0;
 max_adc=20000.;
//=== AC Threshold variable ===//
 ac1_adc[0]=max_ac1;
 ac1_adc[1]=500.;
 ac1_adc[2]=300.;
 ac2_adc[0]=max_ac2;
 ac2_adc[1]=6000.;
 ac2_adc[2]=5000.;
 

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
 ac1_adc[0]=max_ac1;
 ac1_adc[1]=1.3;
 ac1_adc[2]=1.0;
 ac2_adc[0]=max_ac2;
 ac2_adc[1]=13;
 ac2_adc[2]=5;

}


 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 int bin_beta=6000;
 int bin_adc=max_adc-min_adc;
 int bin_ac1=(max_ac1-min_ac1)*3; 
 int bin_ac2=(max_ac2-min_ac2)*3; 

 TH2F* hcoin_ac1[nth];
 TH2F* hcoin_ac2[nth];
 TH1F* hcoin_t1[nth];
 TH1F* hcoin_t2[nth];
 TH1F* hcoin_t3[nth][nth];
 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);

 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);
 
 for(int i=0;i<nth;i++){
 hcoin_ac1[i]=new TH2F(Form("hcoin_ac1[%d]",i),"Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 hcoin_ac2[i]=new TH2F(Form("hcoin_ac2[%d]",i),"Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);

 hcoin_t1[i]=new TH1F(Form("hcoin_t1[%d]",i), Form("Coincidence AC1<%lf cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
 hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence AC2<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);

 for(int j=0;j<nth;j++)hcoin_t3[i][j]=new TH1F(Form("hcoin_t3[%d][%d]",i,j),Form("Coincidence time S2R-S2L[ns] ac1_adc< %lf, ac2_adc<%lf; Coin-Time [ns];Counts ",ac1_adc[i],ac2_adc[j]),bin_coin_c,min_coin_c,max_coin_c);


 } 


 
 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 // double mh;
 double m2; 
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
 double tof_r,tof_l; 
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr,pathl_off;
 pathl_off=0.0;




 // double mm;
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

 if(mode=="G"){
   coin_t=ctime[0];
   coin_tc=ctime[0];
 }else{
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l; //coin time
 coin_tc=coin_t+rpath_corr-lpath_corr+pathl_off;
 
}
//====== Cut condition ========================// 
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_track=false;
   cut_s0=false;
   coin_trig=false;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
  
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


 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){
 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2

 //--------- with AC1 Cut ---------------// 
  for(int j=0;j<nth;j++){
   cut_ac1=false;
   if(Ra1sum<ac1_adc[j])cut_ac1=true;
   if(cut_ac1){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
     }
    }
  //-------with AC2 Cut --------------------//
   for(int k=0;k<nth;k++){
    cut_ac2=false;
    if(Ra2sum<ac2_adc[k])cut_ac2=true;
    if(cut_ac2){
    hcoin_t1[k]->Fill(coin_tc); // AC2 cut
    hcoin_ac1[k]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  
    }
   }
 
 }
 //-------with AC1 && AC2 Cut --------------------//
  cut_ac1=false;
  cut_ac2=false; 
   for(int k=0;k<nth;k++){
   for(int j=0;j<nth;j++){ 
  if(Ra1sum<ac1_adc[j])cut_ac1=true;
  if(Ra2sum<ac2_adc[k])cut_ac2=true;
    
   if(cut_ac1 && cut_ac2)hcoin_t3[j][k]->Fill(coin_tc); //AC1 & AC2 cut
  } 
 }
 
 }

 


 //============== Efficiency analysis =======================//
 
 int iter_ac1=50; int iter_ac2=50; int iter_max=iter_ac1+iter_ac2;
 TF1* facc[iter_max][nth][2];
 TF1* fpi[iter_max][nth][2];
 TF1* fk[iter_max][nth][2];
 TF1* fcoin[iter_max][nth][2];
 TF1* fp[iter_max][nth][2];
 TH1D* hcoin_ac1_p[iter_max][nth];
 TH1D* hcoin_ac2_p[iter_max][nth]; 
 TGraphErrors* gsum_pi_ac1[nth][nth];
 TGraphErrors* gsum_p_ac1[nth][nth];
 TGraphErrors* gsum_k_ac1[nth][nth];
 TGraphErrors* grate_k_ac1[nth][nth];
 TGraphErrors* grate_p_ac1[nth][nth];
 TGraphErrors* grate_pi_ac1[nth][nth];
 TGraphErrors* gsum_pi_ac2[nth][nth];
 TGraphErrors* gsum_p_ac2[nth][nth];
 TGraphErrors* gsum_k_ac2[nth][nth];
 TGraphErrors* grate_k_ac2[nth][nth];
 TGraphErrors* grate_pi_ac2[nth][nth];
 TGraphErrors* grate_p_ac2[nth][nth];
 for(int k=0;k<nth;k++){
 for(int j=0;j<nth;j++){
 gsum_pi_ac1[j][k]=new TGraphErrors();
 gsum_p_ac1[j][k]=new TGraphErrors();
 gsum_k_ac1[j][k]=new TGraphErrors();
 grate_k_ac1[j][k]=new TGraphErrors();
 grate_p_ac1[j][k]=new TGraphErrors();
 grate_pi_ac1[j][k]=new TGraphErrors();
 gsum_pi_ac2[j][k]=new TGraphErrors();
 gsum_p_ac2[j][k]=new TGraphErrors();
 gsum_k_ac2[j][k]=new TGraphErrors();
 grate_k_ac2[j][k]=new TGraphErrors();
 grate_pi_ac2[j][k]=new TGraphErrors();
 grate_p_ac2[j][k]=new TGraphErrors();
 //--- TGraphErrors Setting ----//
 gsum_pi_ac1[j][k]->SetTitle(Form("Sum of Pion (AC2<%lf);AC1 ADC-th [ch];Pion Events [Counts]",ac2_adc[k]));
 gsum_p_ac1[j][k]->SetTitle(Form("SUM of Proton (AC2<%lf);AC1 ADC-th [ch];Proton Events [Counts]",ac2_adc[k]));
 gsum_k_ac1[j][k]->SetTitle(Form("SUM of Kaon (AC2<%lf);AC1 ADC-th [ch];Kaon Events [Counts]",ac2_adc[k]));
 grate_k_ac1[j][k]->SetTitle(Form("Kaoin Survival rate ((AC2<%lf));AC1 ADC-th [ch]; Survival rate ",ac2_adc[k]));
 grate_pi_ac1[j][k]->SetTitle(Form("Pion Survuval rate (AC<%lf);AC1 ADC-th [ch]; Surival rate ",ac2_adc[k]));
 grate_p_ac1[j][k]->SetTitle("Proton Survival rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Survival rate "); 
 gsum_pi_ac1[j][k]->SetMarkerStyle(21);
 gsum_pi_ac1[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac1[j][k]->SetMarkerSize(0.5);
 gsum_p_ac1[j][k]->SetMarkerStyle(21);
 gsum_p_ac1[j][k]->SetMarkerColor(kRed);
 gsum_p_ac1[j][k]->SetMarkerSize(0.5);
 gsum_k_ac1[j][k]->SetMarkerStyle(21);
 gsum_k_ac1[j][k]->SetMarkerColor(kRed);
 gsum_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_k_ac1[j][k]->SetMarkerStyle(21);
 grate_k_ac1[j][k]->SetMarkerColor(kBlue);
 grate_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_p_ac1[j][k]->SetMarkerStyle(21);
 grate_p_ac1[j][k]->SetMarkerColor(kBlue);
 grate_p_ac1[j][k]->SetMarkerSize(0.5);
 grate_pi_ac1[j][k]->SetMarkerStyle(21);
 grate_pi_ac1[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac1[j][k]->SetMarkerSize(0.5);

 gsum_pi_ac2[j][k]->SetTitle("SUM of Pion vs AC2 Threshold;AC2 ADC-th [ch];Pion Events [Counts]");
 gsum_p_ac2[j][k]->SetTitle("SUM of Proton vs AC2 Threshold;AC2 ADC-th [ch];Proton Events [Counts]");
 gsum_k_ac2[j][k]->SetTitle("SUM of Kaon vs AC2 Threshold;AC2 ADC-th [ch];Kaon Events [Counts]");
 grate_k_ac2[j][k]->SetTitle("Kaoin Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_pi_ac2[j][k]->SetTitle("Pion Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_p_ac2[j][k]->SetTitle("Proton Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 gsum_pi_ac2[j][k]->SetMarkerStyle(21);
 gsum_pi_ac2[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac2[j][k]->SetMarkerSize(0.5);
 gsum_p_ac2[j][k]->SetMarkerStyle(21);
 gsum_p_ac2[j][k]->SetMarkerColor(kRed);
 gsum_p_ac2[j][k]->SetMarkerSize(0.5);
 gsum_k_ac2[j][k]->SetMarkerStyle(21);
 gsum_k_ac2[j][k]->SetMarkerColor(kRed);
 gsum_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_k_ac2[j][k]->SetMarkerStyle(21);
 grate_k_ac2[j][k]->SetMarkerColor(kBlue);
 grate_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_p_ac2[j][k]->SetMarkerStyle(21);
 grate_p_ac2[j][k]->SetMarkerColor(kBlue);
 grate_p_ac2[j][k]->SetMarkerSize(0.5);
 grate_pi_ac2[j][k]->SetMarkerStyle(21);
 grate_pi_ac2[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac2[j][k]->SetMarkerSize(0.5);
 }
 }
//--- Parameters -----//

 double kmin[iter_max][nth][2],kmax[iter_max][nth][2];
 double inte_ktot[iter_max][nth][2], inte_ksig[iter_max][nth][2];
 double p0_acc[iter_max][nth][2], p1_acc[iter_max][nth][2];
 double n_p[iter_max][nth][2],sig_p[iter_max][nth][2],mean_p[iter_max][nth][2];
 double n_pi[iter_max][nth][2],sig_pi[iter_max][nth][2],mean_pi[iter_max][nth][2];
 double n_k[iter_max][nth][2],sig_k[iter_max][nth][2],mean_k[iter_max][nth][2];
 int bin_ac1_adc[nth][nth],bin_min_ac1,bin_ac2_adc[nth][nth],bin_max_ac2,bin_min_ac2;
 double sum_k[iter_max][nth][2],sum_p[iter_max][nth][2],sum_pi[iter_max][nth][2]; 
  double sum_k_err[iter_max][nth][2],sum_p_err[iter_max][nth][2],sum_pi_err[iter_max][nth][2]; 
double inte_acc[iter_max][nth][2];
 double th_ac1[iter_max],th_ac2[iter_max];
 int bin_th_ac1[iter_max][nth],bin_th_ac2[iter_max][nth]; 
 double nk[iter_max][nth][iter_max][nth][2],npi[iter_max][nth][iter_max][nth][2],np[iter_max][nth][iter_max][nth][2];
 double max_nk[nth][nth][2],max_npi[nth][nth][2],max_np[nth][nth][2];
 double n_p_err[iter_max][nth][2],n_pi_err[iter_max][nth][2],n_k_err[iter_max][nth][2];
//---- Defolt parameters -----------//

 TF1* facc_t1def[nth][nth];
 TF1* fpi_t1def[nth][nth];
 TF1* fk_t1def[nth][nth];
 TF1* fcoin_t1def[nth][nth];
 TF1* fp_t1def[nth][nth];
 TF1* facc_t2def[nth][nth];
 TF1* fpi_t2def[nth][nth];
 TF1* fk_t2def[nth][nth];
 TF1* fcoin_t2def[nth][nth];
 TF1* fp_t2def[nth][nth];
 TF1* facc_t3def[nth][nth];
 TF1* fpi_t3def[nth][nth];
 TF1* fk_t3def[nth][nth];
 TF1* fcoin_t3def[nth][nth];
 TF1* fp_t3def[nth][nth];
 TF1* fcoin_t1[nth][nth];
 TF1* fcoin_t2[nth][nth];  
 TF1* fcoin_t3[nth][nth]; 


 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_sig_k,def_mean_k,def_acc;


 if(mode=="H"){
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


 double def_t1_k[nth][nth],def_t1_pi[nth][nth],def_t1_p[nth][nth],def_t1_acc[nth][nth];
 double def_t1_k_err[nth][nth],def_t1_pi_err[nth][nth],def_t1_p_err[nth][nth],def_t1_acc_err[nth][nth];
 double t1sig_k[nth][nth],t1sig_p[nth][nth],t1sig_pi[nth][nth],t1mean_p[nth][nth],t1mean_k[nth][nth],t1mean_pi[nth][nth];
 double t1sum_k[nth],t1sum_pi[nth],t1sum_p[nth];
double t1sum_k_err[nth],t1sum_pi_err[nth],t1sum_p_err[nth];
 double def_t2_k[nth][nth],def_t2_pi[nth][nth],def_t2_p[nth][nth],def_t2_acc[nth][nth];
 double t2sig_k[nth][nth],t2sig_p[nth][nth],t2sig_pi[nth][nth],t2mean_p[nth][nth],t2mean_k[nth][nth],t2mean_pi[nth][nth];
 double def_t2_k_err[nth][nth],def_t2_pi_err[nth][nth],def_t2_p_err[nth][nth],def_t2_acc_err[nth][nth];
 double t2sum_k[nth],t2sum_pi[nth],t2sum_p[nth];
 double t2sum_k_err[nth],t2sum_pi_err[nth],t2sum_p_err[nth];
 double def_t3_k[nth][nth],def_t3_pi[nth][nth],def_t3_p[nth][nth],def_t3_acc[nth][nth];
 double t3sig_k[nth][nth],t3sig_p[nth][nth],t3sig_pi[nth][nth],t3mean_p[nth][nth],t3mean_k[nth][nth],t3mean_pi[nth][nth];
 double t3sum_k[nth][nth],t3sum_pi[nth][nth],t3sum_p[nth][nth];
 double emp[iter_max];
 for(int i;i<iter_max;i++){emp[i]=0.0;}
 bool point_err=true;
//----------------------------------//
 



      for(int th1=0;th1<nth;th1++){
	for(int th2=0;th2<nth;th2++){



	 

   //========================================================//
  //========== Get Maximum Events ===========================//
   //=======================================================//

//-----Defined Function ---------//

facc_t1def[th1][th2]=new TF1(Form("facc_t1def[%d][%d]",th1,th2),"[0]",min_coin_c,max_coin_c);
fp_t1def[th1][th2] =new TF1(Form("fp_t1def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fpi_t1def[th1][th2]=new TF1(Form("fpi_t1def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fk_t1def[th1][th2] =new TF1(Form("fk_t1def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fcoin_t1[th1][th2]=new TF1(Form("fcoin_t1[%d][%d]",th1,th2),"gausn(0)+gausn(3)+gausn(6)+pol1(9)"); 
facc_t2def[th1][th2]=new TF1(Form("facc_t2def[%d][%d]",th1,th2),"[0]",min_coin_c,max_coin_c);
fp_t2def[th1][th2] =new TF1(Form("fp_t2def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fpi_t2def[th1][th2]=new TF1(Form("fpi_t2def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fk_t2def[th1][th2] =new TF1(Form("fk_t2def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fcoin_t2[th1][th2]=new TF1(Form("fcoin_t2[%d][%d]",th1,th2),"gausn(0)+gausn(3)+gausn(6)+pol1(9)"); 
facc_t3def[th1][th2]=new TF1(Form("facc_t3def[%d][%d]",th1,th2),"[0]",min_coin_c,max_coin_c);
fp_t3def[th1][th2] =new TF1(Form("fp_t3def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fpi_t3def[th1][th2]=new TF1(Form("fpi_t3def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fk_t3def[th1][th2] =new TF1(Form("fk_t3def[%d][%d]",th1,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
fcoin_t3[th1][th2]=new TF1(Form("fcoin_t3[%d][%d]",th1,th2),"gausn(0)+gausn(3)+gausn(6)+pol1(9)"); 
//--------- Fitting Function ------//



 hcoin_t1[th2]->Fit(Form("facc_t1def[%d][%d]",th1,th2),"R","",min_coin_c,def_mean_p-5*def_sig_p);
 def_t1_acc[th1][th2]=facc_t1def[th1][th2]->GetParameter(0);
 fpi_t1def[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi_t1def[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi_t1def[th1][th2]->FixParameter(3,def_t1_acc[th1][th2]);
 fk_t1def[th1][th2]->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk_t1def[th1][th2]->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk_t1def[th1][th2]->FixParameter(3,def_t1_acc[th1][th2]);
 fp_t1def[th1][th2]->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp_t1def[th1][th2]->SetParLimits(2,0.5*def_sig_p,1.5*def_sig_p);
 fp_t1def[th1][th2]->FixParameter(3,def_t1_acc[th1][th2]);
 hcoin_t1[th2]->Fit(Form("facc_t1def[%d][%d]",th1,th2),"Rq","",min_coin_c,def_mean_p-5*def_sig_p);
 hcoin_t1[th2]->Fit(Form("fp_t1def[%d][%d]",th1,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 hcoin_t1[th2]->Fit(Form("fpi_t1def[%d][%d]",th1,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_t1[th2]->Fit(Form("fk_t1def[%d][%d]",th1,th2),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);


 hcoin_t2[th1]->Fit(Form("facc_t2def[%d][%d]",th1,th2),"Rq","",min_coin_c,def_mean_p-5*def_sig_p);
 def_t2_acc[th1][th2]=facc_t2def[th1][th2]->GetParameter(0);
 fpi_t2def[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi_t2def[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi_t2def[th1][th2]->FixParameter(3,def_t2_acc[th1][th2]);
 fk_t2def[th1][th2]->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk_t2def[th1][th2]->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk_t2def[th1][th2]->FixParameter(3,def_t2_acc[th1][th2]);
 fp_t2def[th1][th2]->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp_t2def[th1][th2]->SetParLimits(2,0.5*def_sig_p,1.5*def_sig_p);
 fp_t2def[th1][th2]->FixParameter(3,def_t2_acc[th1][th2]);
 hcoin_t2[th1]->Fit(Form("facc_t2def[%d][%d]",th1,th2),"Rq","",min_coin_c,def_mean_p-5*def_sig_p);
 hcoin_t2[th1]->Fit(Form("fp_t2def[%d][%d]",th1,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 hcoin_t2[th1]->Fit(Form("fpi_t2def[%d][%d]",th1,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_t2[th1]->Fit(Form("fk_t2def[%d][%d]",th1,th2),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);

 hcoin_t3[th1][th2]->Fit(Form("facc_t3def[%d][%d]",th1,th2),"Rq","",min_coin_c,def_mean_p-5*def_sig_p);
 def_t3_acc[th1][th2]=facc_t3def[th1][th2]->GetParameter(0);
 fpi_t3def[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi_t3def[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi_t3def[th1][th2]->FixParameter(3,def_t3_acc[th1][th2]);
 fk_t3def[th1][th2]->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk_t3def[th1][th2]->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk_t3def[th1][th2]->FixParameter(3,def_t3_acc[th1][th2]);
 fp_t3def[th1][th2]->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp_t3def[th1][th2]->SetParLimits(2,0.5*def_sig_p,1.5*def_sig_p);
 fp_t3def[th1][th2]->FixParameter(3,def_t3_acc[th1][th2]);
 hcoin_t3[th1][th2]->Fit(Form("facc_t3def[%d][%d]",th1,th2),"Rq","",min_coin_c,def_mean_p-5*def_sig_p);
 hcoin_t3[th1][th2]->Fit(Form("fp_t3def[%d][%d]",th1,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 hcoin_t3[th1][th2]->Fit(Form("fpi_t3def[%d][%d]",th1,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_t3[th1][th2]->Fit(Form("fk_t3def[%d][%d]",th1,th2),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);



 //------- Getting Initial Parameters ----------//

 def_t1_p[th1][th2]=fp_t1def[th1][th2]->GetParameter(0); t1mean_p[th1][th2]=fp_t1def[th1][th2]->GetParameter(1);  t1sig_p[th1][th2]=fp_t1def[th1][th2]->GetParameter(2);  
 def_t1_pi[th1][th2]=fpi_t1def[th1][th2]->GetParameter(0); t1mean_pi[th1][th2]=fpi_t1def[th1][th2]->GetParameter(1);  t1sig_pi[th1][th2]=fpi_t1def[th1][th2]->GetParameter(2); 
 def_t1_k[th1][th2]=fk_t1def[th1][th2]->GetParameter(0); t1mean_k[th1][th2]=fk_t1def[th1][th2]->GetParameter(1);  t1sig_k[th1][th2]=fk_t1def[th1][th2]->GetParameter(2);

 def_t2_p[th1][th2]=fp_t2def[th1][th2]->GetParameter(0); t2mean_p[th1][th2]=fp_t2def[th1][th2]->GetParameter(1);  t2sig_p[th1][th2]=fp_t2def[th1][th2]->GetParameter(2);  
 def_t2_pi[th1][th2]=fpi_t2def[th1][th2]->GetParameter(0); t2mean_pi[th1][th2]=fpi_t2def[th1][th2]->GetParameter(1);  t2sig_pi[th1][th2]=fpi_t2def[th1][th2]->GetParameter(2); 
 def_t2_k[th1][th2]=fk_t2def[th1][th2]->GetParameter(0); t2mean_k[th1][th2]=fk_t2def[th1][th2]->GetParameter(1);  t2sig_k[th1][th2]=fk_t2def[th1][th2]->GetParameter(2);


 def_t3_p[th1][th2]=fp_t3def[th1][th2]->GetParameter(0); t3mean_p[th1][th2]=fp_t3def[th1][th2]->GetParameter(1);  t3sig_p[th1][th2]=fp_t3def[th1][th2]->GetParameter(2);  
 def_t3_pi[th1][th2]=fpi_t3def[th1][th2]->GetParameter(0); t3mean_pi[th1][th2]=fpi_t3def[th1][th2]->GetParameter(1);  t3sig_pi[th1][th2]=fpi_t3def[th1][th2]->GetParameter(2); 
 def_t3_k[th1][th2]=fk_t3def[th1][th2]->GetParameter(0); t3mean_k[th1][th2]=fk_t3def[th1][th2]->GetParameter(1);  t3sig_k[th1][th2]=fk_t3def[th1][th2]->GetParameter(2);


 //--------- coin  Set Paramters -----------//


 fcoin_t1[th1][th2]->SetParameters(def_t1_pi[th1][th2],t1mean_pi[th1][th2],t1sig_pi[th1][th2],def_t1_k[th1][th2],t1mean_k[th1][th2],t1sig_k[th1][th2],def_t1_p[th1][th2],t1mean_p[th1][th2],t1sig_p[th1][th2],def_t1_acc[th1][th2]);// Parameters =Pion(0) ,TH2aon(3) , Proton(6) , acc(9)
 fcoin_t1[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t1[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t1[th1][th2]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t1[th1][th2]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin_t1[th1][th2]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t1[th1][th2]->SetParLimits(8,0.5*def_sig_p,1.5*def_sig_p);
 fcoin_t1[th1][th2]->FixParameter(9,def_t1_acc[th1][th2]);



 fcoin_t2[th1][th2]->SetParameters(def_t2_pi[th1][th2],t2mean_pi[th1][th2],t2sig_pi[th1][th2],def_t2_k[th1][th2],t2mean_k[th1][th2],t2sig_k[th1][th2],def_t2_p[th1][th2],t2mean_p[th1][th2],t2sig_p[th1][th2],def_t2_acc[th1][th2]);// Parameters =Pion(0) ,TH2aon(3) , Proton(6) , acc(9)
 fcoin_t2[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t2[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t2[th1][th2]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t2[th1][th2]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin_t2[th1][th2]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t2[th1][th2]->SetParLimits(8,0.5*def_sig_p,1.5*def_sig_p);
 fcoin_t2[th1][th2]->FixParameter(9,def_t2_acc[th1][th2]);




 fcoin_t3[th1][th2]->SetParameters(def_t3_pi[th1][th2],t3mean_pi[th1][th2],t3sig_pi[th1][th2],def_t3_k[th1][th2],t3mean_k[th1][th2],t3sig_k[th1][th2],def_t3_p[th1][th2],t3mean_p[th1][th2],t3sig_p[th1][th2],def_t3_acc[th1][th2]);// Parameters =Pion(0) ,TH2aon(3) , Proton(6) , acc(9)
 fcoin_t3[th1][th2]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin_t3[th1][th2]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin_t3[th1][th2]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin_t3[th1][th2]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin_t3[th1][th2]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin_t3[th1][th2]->SetParLimits(8,0.5*def_sig_p,1.5*def_sig_p);
 fcoin_t3[th1][th2]->FixParameter(9,def_t3_acc[th1][th2]);



 //------ Fitting Coincidence Function -------//
 hcoin_t1[th2]->Fit(Form("fcoin_t1[%d][%d]",th1,th2),"Rq","",min_coin_c,max_coin_c);
 hcoin_t2[th1]->Fit(Form("fcoin_t2[%d][%d]",th1,th2),"Rq","",min_coin_c,max_coin_c);
 hcoin_t3[th1][th2]->Fit(Form("fcoin_t3[%d][%d]",th1,th2),"Rq","",min_coin_c,max_coin_c);
 //------ Getting Parameters -------------//

 def_t1_pi[th1][th2]=fcoin_t1[th1][th2]->GetParameter(0);t1mean_pi[th1][th2]=fcoin_t1[th1][th2]->GetParameter(1); t1sig_pi[th1][th2]=fcoin_t1[th1][th2]->GetParameter(2);  
 def_t1_k[th1][th2]=fcoin_t1[th1][th2]->GetParameter(3); t1mean_k[th1][th2]=fcoin_t1[th1][th2]->GetParameter(4);  t1sig_k[th1][th2]=fcoin_t1[th1][th2]->GetParameter(5);  
 def_t1_p[th1][th2]=fcoin_t1[th1][th2]->GetParameter(6); t1mean_p[th1][th2]=fcoin_t1[th1][th2]->GetParameter(7);  t1sig_p[th1][th2]=fcoin_t1[th1][th2]->GetParameter(8);  

 def_t2_pi[th1][th2]=fcoin_t2[th1][th2]->GetParameter(0);t2mean_pi[th1][th2]=fcoin_t2[th1][th2]->GetParameter(1); t2sig_pi[th1][th2]=fcoin_t2[th1][th2]->GetParameter(2);  
 def_t2_k[th1][th2]=fcoin_t2[th1][th2]->GetParameter(3); t2mean_k[th1][th2]=fcoin_t2[th1][th2]->GetParameter(4);  t2sig_k[th1][th2]=fcoin_t2[th1][th2]->GetParameter(5);  
 def_t2_p[th1][th2]=fcoin_t2[th1][th2]->GetParameter(6); t2mean_p[th1][th2]=fcoin_t2[th1][th2]->GetParameter(7);  t2sig_p[th1][th2]=fcoin_t2[th1][th2]->GetParameter(8);  

 def_t3_pi[th1][th2]=fcoin_t3[th1][th2]->GetParameter(0);t3mean_pi[th1][th2]=fcoin_t3[th1][th2]->GetParameter(1); t3sig_pi[th1][th2]=fcoin_t3[th1][th2]->GetParameter(2);  
 def_t3_k[th1][th2]=fcoin_t3[th1][th2]->GetParameter(3); t3mean_k[th1][th2]=fcoin_t3[th1][th2]->GetParameter(4);  t3sig_k[th1][th2]=fcoin_t3[th1][th2]->GetParameter(5);  
 def_t3_p[th1][th2]=fcoin_t3[th1][th2]->GetParameter(6); t3mean_p[th1][th2]=fcoin_t3[th1][th2]->GetParameter(7);  t3sig_p[th1][th2]=fcoin_t3[th1][th2]->GetParameter(8);  

 //------ Getting Error Paramters --------//
 def_t1_pi_err[th1][th2]=fcoin_t1[th1][th2]->GetParError(0);
 def_t1_k_err[th1][th2]=fcoin_t1[th1][th2]->GetParError(3); 
 def_t1_p_err[th1][th2]=fcoin_t1[th1][th2]->GetParError(6); 
 def_t2_pi_err[th1][th2]=fcoin_t2[th1][th2]->GetParError(0);
 def_t2_k_err[th1][th2]=fcoin_t2[th1][th2]->GetParError(3);   
 def_t2_p_err[th1][th2]=fcoin_t2[th1][th2]->GetParError(6); 


 //---- Set Paramters ------------//

 fpi_t1def[th1][th2]->SetParameter(0,def_t1_pi[th1][th2]);
 fpi_t1def[th1][th2]->SetParameter(1,t1mean_pi[th1][th2]);
 fpi_t1def[th1][th2]->SetParameter(2,t1sig_pi[th1][th2]);
 fpi_t1def[th1][th2]->SetParameter(3,def_t1_acc[th1][th2]);
 fp_t1def[th1][th2]->SetParameter(0,def_t1_p[th1][th2]);
 fp_t1def[th1][th2]->SetParameter(1,t1mean_p[th1][th2]);
 fp_t1def[th1][th2]->SetParameter(2,t1sig_p[th1][th2]);
 fp_t1def[th1][th2]->SetParameter(3,def_t1_acc[th1][th2]);
 fk_t1def[th1][th2]->SetParameter(0,def_t1_k[th1][th2]);
 fk_t1def[th1][th2]->SetParameter(1,t1mean_k[th1][th2]);
 fk_t1def[th1][th2]->SetParameter(2,t1sig_k[th1][th2]);
 fk_t1def[th1][th2]->SetParameter(3,def_t1_acc[th1][th2]);

 fpi_t2def[th1][th2]->SetParameter(0,def_t2_pi[th1][th2]);
 fpi_t2def[th1][th2]->SetParameter(1,t2mean_pi[th1][th2]);
 fpi_t2def[th1][th2]->SetParameter(2,t2sig_pi[th1][th2]);
 fpi_t2def[th1][th2]->SetParameter(3,def_t2_acc[th1][th2]);
 fp_t2def[th1][th2]->SetParameter(0,def_t2_p[th1][th2]);
 fp_t2def[th1][th2]->SetParameter(1,t2mean_p[th1][th2]);
 fp_t2def[th1][th2]->SetParameter(2,t2sig_p[th1][th2]);
 fp_t2def[th1][th2]->SetParameter(3,def_t2_acc[th1][th2]);
 fk_t2def[th1][th2]->SetParameter(0,def_t2_k[th1][th2]);
 fk_t2def[th1][th2]->SetParameter(1,t2mean_k[th1][th2]);
 fk_t2def[th1][th2]->SetParameter(2,t2sig_k[th1][th2]);
 fk_t2def[th1][th2]->SetParameter(3,def_t2_acc[th1][th2]);

 fpi_t3def[th1][th2]->SetParameter(0,def_t3_pi[th1][th2]);
 fpi_t3def[th1][th2]->SetParameter(1,t3mean_pi[th1][th2]);
 fpi_t3def[th1][th2]->SetParameter(2,t3sig_pi[th1][th2]);
 fpi_t3def[th1][th2]->SetParameter(3,def_t3_acc[th1][th2]);
 fp_t3def[th1][th2]->SetParameter(0,def_t3_p[th1][th2]);
 fp_t3def[th1][th2]->SetParameter(1,t3mean_p[th1][th2]);
 fp_t3def[th1][th2]->SetParameter(2,t3sig_p[th1][th2]);
 fp_t3def[th1][th2]->SetParameter(3,def_t3_acc[th1][th2]);
 fk_t3def[th1][th2]->SetParameter(0,def_t3_k[th1][th2]);
 fk_t3def[th1][th2]->SetParameter(1,t3mean_k[th1][th2]);
 fk_t3def[th1][th2]->SetParameter(2,t3sig_k[th1][th2]);
 fk_t3def[th1][th2]->SetParameter(3,def_t3_acc[th1][th2]);



 //======== SUM[ac1_adc] ==========//
t1sum_k[th2]=def_t1_k[th1][th2]/tdc_time;
t1sum_p[th2]=def_t1_p[th1][th2]/tdc_time;
t1sum_pi[th2]=def_t1_pi[th1][th2]/tdc_time;

t1sum_k_err[th2]=def_t1_k_err[th1][th2]/tdc_time;
t1sum_p_err[th2]=def_t1_p_err[th1][th2]/tdc_time;
t1sum_pi_err[th2]=def_t1_pi_err[th1][th2]/tdc_time;
 //======== SUM[ac2_adc] ==========//
t2sum_k[th1]=def_t2_k[th1][th2]/tdc_time;
t2sum_p[th1]=def_t2_p[th1][th2]/tdc_time;
t2sum_pi[th1]=def_t2_pi[th1][th2]/tdc_time;

t2sum_k_err[th1]=def_t2_k_err[th1][th2]/tdc_time;
t2sum_p_err[th1]=def_t2_p_err[th1][th2]/tdc_time;
t2sum_pi_err[th1]=def_t2_pi_err[th1][th2]/tdc_time;
//===== SUM[ac1_adc][ac2_adc] ==========//
t3sum_k[th1][th2]=def_t3_k[th1][th2]/tdc_time;
t3sum_p[th1][th2]=def_t3_p[th1][th2]/tdc_time;
t3sum_pi[th1][th2]=def_t3_pi[th1][th2]/tdc_time;

	 



 max_nk[th1][th2][0]=0.0; max_nk[th1][th2][1]=0.0; max_np[th1][th2][0]=0.0; max_np[th1][th2][1]=0.0; 
 max_npi[th1][th2][0]=0.0; max_npi[th1][th2][1]=0.0;


 //=======================================================================================//




     for(int i=0;i<iter_ac1;i++){

       
       // th_ac1[i]=min_ac1+(max_ac1-min_ac1)/iter_ac1*i;
       // th_ac2[i]=min_ac2+(max_ac2-min_ac2)/iter_ac1*i;	

          th_ac1[i]=max_ac1-(max_ac1-min_ac1)/iter_ac1*i;
          th_ac2[i]=max_ac2-(max_ac2-min_ac2)/iter_ac1*i;	

 bin_th_ac1[i][th2]=hcoin_ac1[th2]->GetYaxis()->FindBin(th_ac1[i]);//hcoin_ac1[th2] is ac2 fixed threshold
 bin_min_ac1=hcoin_ac1[th2]->GetYaxis()->FindBin(min_ac1);
 hcoin_ac1_p[i][th2]=hcoin_ac1[th2]->ProjectionX(Form("hcoin_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);

 bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
 bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
 hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); //hcoin_ac2[th1] is ac1 fixed threshold

 //--- Initial Parameters -----------//

 facc[i][th2][0]=new TF1(Form("facc[%d][%d][0]",i,th2),"[0]",min_coin_c,max_coin_c);
  hcoin_ac1_p[i][th2]->Fit(Form("facc[%d][%d][0]",i,th2),"Rq","",min_coin_c,min_coin_c+3.);
 p0_acc[i][th2][0]=facc[i][th2][0]->GetParameter(0);
 
facc[i][th1][1]=new TF1(Form("facc[%d][%d][1]",i,th1),"[0]",min_coin_c,max_coin_c);  
hcoin_ac2_p[i][th1]->Fit(Form("facc[%d][%d][1]",i,th1),"Rq","",min_coin_c,min_coin_c+3.);
 p0_acc[i][th1][1]=facc[i][th1][1]->GetParameter(0);

 //------- AC1 --------------// 

 fp[i][th2][0]=new TF1(Form("fp[%d][%d][0]",i,th2),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi[i][th2][0] =new TF1(Form("fpi[%d][%d][0]",i,th2),"gausn+pol0(3)",min_coin_c,max_coin_c);
 fk[i][th2][0]=new TF1(Form("fk[%d][%d][0]",i,th2),"gausn+pol0(3)",min_coin_c,max_coin_c);
 fp[i][th2][0]->SetParameter(1,def_mean_p); fp[i][th2][0]->SetParameter(2,def_sig_p); fp[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);
 fpi[i][th2][0]->SetParameter(1,def_mean_pi); fpi[i][th2][0]->SetParameter(2,def_sig_pi); fpi[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);
 fk[i][th2][0]->SetParameter(1,def_mean_k); fk[i][th2][0]->SetParameter(2,def_sig_k);fk[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);
 fcoin[i][th2][0] =new TF1(Form("fcoin[%d][%d][0]",i,th2),"gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin[i][th2][0]->SetTitle(Form("Coin-Time w AC cut  (AC1<%d ch && AC2>%d ch);Coin time [ns];Counts [1/56 ns]",th_ac1[i],ac2_adc));
 //------- AC2 -------------//

 fp[i][th1][1]=new TF1(Form("fp[%d][%d][1]",i,th1),"gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi[i][th1][1] =new TF1(Form("fpi[%d][%d][1]",i,th1),"gausn+pol0(3)",min_coin_c,max_coin_c);
 fk[i][th1][1]=new TF1(Form("fk[%d][%d][1]",i,th1),"gausn+pol0(3)",min_coin_c,max_coin_c);
 fpi[i][th1][1]->SetParameter(1,def_mean_pi); fpi[i][th1][1]->SetParameter(2,def_sig_pi); fpi[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);
 fk[i][th1][1]->SetParameter(1,def_mean_k); fk[i][th1][1]->SetParameter(2,def_sig_k); fk[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);
 fp[i][th1][1]->SetParameter(1,def_mean_p); fp[i][th1][1]->SetParameter(2,def_sig_p); fp[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);

 fcoin[i][th1][1] =new TF1(Form("fcoin[%d][%d][1]",i,th1),"gausn(0)+gausn(3)+gausn(6)+pol1(9)");
 fcoin[i][th1][1]->SetTitle(Form("Coin-Time w AC cut (AC1 Cut<%lf ch && AC2 Cut>%lf ch);Coin time [ns];Counts [1/56 ns]",ac1_adc,th_ac2[i]));
 //----------------------------------//

 
 //----- AC1 Fitting -----------// 
 hcoin_ac1_p[i][th2]->Fit(Form("fp[%d][%d][0]",i,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th2][0]=fp[i][th2][0]->GetParameter(0);
 mean_p[i][th2][0]=fp[i][th2][0]->GetParameter(1);
 sig_p[i][th2][0]=fp[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fpi[%d][%d][0]",i,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(0);
 mean_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(1);
 sig_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fk[%d][%d][0]",i,th2),"Rq","",def_mean_p-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th2][0]=fk[i][th2][0]->GetParameter(0);
 mean_k[i][th2][0]=fk[i][th2][0]->GetParameter(1);
 sig_k[i][th2][0]=fk[i][th2][0]->GetParameter(2);
  //----- AC2 Fitting -----------// 
 hcoin_ac2_p[i][th1]->Fit(Form("fp[%d][%d][1]",i,th1),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th1][1]=fp[i][th1][1]->GetParameter(0);
 mean_p[i][th1][1]=fp[i][th1][1]->GetParameter(1);
 sig_p[i][th1][1]=fp[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fpi[%d][%d][1]",i,th1),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(0);
 mean_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(1);
 sig_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fk[%d][%d][1]",i,th1),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th1][1]=fk[i][th1][1]->GetParameter(0);
 mean_k[i][th1][1]=fk[i][th1][1]->GetParameter(1);
 sig_k[i][th1][1]=fk[i][th1][1]->GetParameter(2);

 //----- AC1 Coint Fitting ---------//

 fcoin[i][th2][0]->SetParameters(n_pi[i][th2][0],mean_pi[i][th2][0],sig_pi[i][th2][0],n_k[i][th2][0],mean_k[i][th2][0],sig_k[i][th2][0],n_p[i][th2][0],mean_p[i][th2][0],sig_p[i][th2][0],p0_acc[i][th2][0]);
 fcoin[i][th2][0]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin[i][th2][0]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin[i][th2][0]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin[i][th2][0]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin[i][th2][0]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin[i][th2][0]->SetParLimits(8,0.8*def_sig_p,1.2*def_sig_p);
 hcoin_ac1_p[i][th2]->Fit(Form("fcoin[%d][%d][0]",i,th2),"Rq","",min_coin_c,max_coin_c);
 
 n_pi[i][th2][0]=fcoin[i][th2][0]->GetParameter(0); 
 mean_pi[i][th2][0]=fcoin[i][th2][0]->GetParameter(1);  
 sig_pi[i][th2][0]=fcoin[i][th2][0]->GetParameter(2);
 n_k[i][th2][0]=fcoin[i][th2][0]->GetParameter(3); 
 mean_k[i][th2][0]=fcoin[i][th2][0]->GetParameter(4);  
 sig_k[i][th2][0]=fcoin[i][th2][0]->GetParameter(5);
 n_p[i][th2][0]=fcoin[i][th2][0]->GetParameter(6); 
 mean_p[i][th2][0]=fcoin[i][th2][0]->GetParameter(7);  
 sig_p[i][th2][0]=fcoin[i][th2][0]->GetParameter(8); 
 p0_acc[i][th2][0]=fcoin[i][th2][0]->GetParameter(9);

 //------- Get Error Paramters ---//

 n_pi_err[i][th2][0]=fcoin[i][th2][0]->GetParError(0); 
 n_k_err[i][th2][0]=fcoin[i][th2][0]->GetParError(3); 
 n_p_err[i][th2][0]=fcoin[i][th2][0]->GetParError(6); 

fp[i][th2][0]->SetParameter(0,n_p[i][th2][0]); fp[i][th2][0]->SetParameter(1,mean_p[i][th2][0]);fp[i][th2][0]->SetParameter(2,sig_p[i][th2][0]); fp[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);
fpi[i][th2][0]->SetParameter(0,n_pi[i][th2][0]); fpi[i][th2][0]->SetParameter(1,mean_pi[i][th2][0]);fpi[i][th2][0]->SetParameter(2,sig_pi[i][th2][0]); fpi[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);
fk[i][th2][0]->SetParameter(0,n_k[i][th2][0]); fk[i][th2][0]->SetParameter(1,mean_k[i][th2][0]);fk[i][th2][0]->SetParameter(2,sig_k[i][th2][0]); fk[i][th2][0]->SetParameter(3,p0_acc[i][th2][0]);


 //----- AC2 Coint Fitting ---------//

 fcoin[i][th1][1]->SetParameters(n_pi[i][th1][1],mean_pi[i][th1][1],sig_pi[i][th1][1],n_k[i][th1][1],mean_k[i][th1][1],sig_k[i][th1][1],n_p[i][th1][1],mean_p[i][th1][1],sig_p[i][th1][1],p0_acc[i][th1][1]);
 fcoin[i][th1][1]->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fcoin[i][th1][1]->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fcoin[i][th1][1]->SetParLimits(4,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fcoin[i][th1][1]->SetParLimits(5,0.8*def_sig_k,1.2*def_sig_k);
 fcoin[i][th1][1]->SetParLimits(7,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fcoin[i][th1][1]->SetParLimits(8,0.8*def_sig_p,1.2*def_sig_p);
 hcoin_ac2_p[i][th1]->Fit(Form("fcoin[%d][%d][1]",i,th1),"Rq","",min_coin_c,max_coin_c);


 n_pi[i][th1][1]=fcoin[i][th1][1]->GetParameter(0); 
 mean_pi[i][th1][1]=fcoin[i][th1][1]->GetParameter(1);  
 sig_pi[i][th1][1]=fcoin[i][th1][1]->GetParameter(2);
 n_k[i][th1][1]=fcoin[i][th1][1]->GetParameter(3); 
 mean_k[i][th1][1]=fcoin[i][th1][1]->GetParameter(4);  
 sig_k[i][th1][1]=fcoin[i][th1][1]->GetParameter(5);
 n_p[i][th1][1]=fcoin[i][th1][1]->GetParameter(6); 
 mean_p[i][th1][1]=fcoin[i][th1][1]->GetParameter(7);  
 sig_p[i][th1][1]=fcoin[i][th1][1]->GetParameter(8);
 p0_acc[i][th1][1]=fcoin[i][th1][1]->GetParameter(9);


 n_pi_err[i][th2][1]=fcoin[i][th1][1]->GetParError(0); 
 n_k_err[i][th2][1]=fcoin[i][th1][1]->GetParError(3); 
 n_p_err[i][th2][1]=fcoin[i][th1][1]->GetParError(6); 


fp[i][th1][1]->SetParameter(0,n_p[i][th1][1]); 
fp[i][th1][1]->SetParameter(1,mean_p[i][th1][1]);
fp[i][th1][1]->SetParameter(2,sig_p[i][th1][1]); 
fp[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);
fpi[i][th1][1]->SetParameter(0,n_pi[i][th1][1]); 
fpi[i][th1][1]->SetParameter(1,mean_pi[i][th1][1]);
fpi[i][th1][1]->SetParameter(2,sig_pi[i][th1][1]); 
fpi[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);
fk[i][th1][1]->SetParameter(0,n_k[i][th1][1]); 
fk[i][th1][1]->SetParameter(1,mean_k[i][th1][1]);
fk[i][th1][1]->SetParameter(2,sig_k[i][th1][1]); 
fk[i][th1][1]->SetParameter(3,p0_acc[i][th1][1]);

  //---- AC1 ----//
 sum_k[i][th2][0]=n_k[i][th2][0]/tdc_time; 
 sum_pi[i][th2][0]=n_pi[i][th2][0]/tdc_time;
 sum_p[i][th2][0]=n_p[i][th2][0]/tdc_time;

 sum_k_err[i][th2][0]=n_k_err[i][th2][0]/tdc_time; 
 sum_pi_err[i][th2][0]=n_pi_err[i][th2][0]/tdc_time;
 sum_p_err[i][th2][0]=n_p_err[i][th2][0]/tdc_time;
 //---- AC2 ----//
 sum_k[i][th1][1]=n_k[i][th1][1]/tdc_time; 
 sum_pi[i][th1][1]=n_pi[i][th1][1]/tdc_time;
 sum_p[i][th1][1]=n_p[i][th1][1]/tdc_time;

 sum_k_err[i][th1][1]=n_k_err[i][th1][1]/tdc_time; 
 sum_pi_err[i][th1][1]=n_pi_err[i][th1][1]/tdc_time;
 sum_p_err[i][th1][1]=n_p_err[i][th1][1]/tdc_time;
 

 //---- AC1 ----//
 gsum_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_pi[i][th2][0]);
 gsum_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_p[i][th2][0]);
 gsum_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]);

 //---- AC2 ----//
 gsum_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_pi[i][th1][1]);
 gsum_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_p[i][th1][1]);
 gsum_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]);



 /*
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]/t1sum_k[th2]);
 grate_k_ac1[th1][th2]->SetPointError(i,emp[i],sum_k[i][th2][0]/t1sum_k[th2]*(t1sum_k_err[th2]/t1sum_k[th2]+sum_k_err[i][th2][0]/sum_k[i][th2][0]));
 grate_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_p[i][th2][0]/t1sum_p[th2]);
 grate_p_ac1[th1][th2]->SetPointError(i,emp[i],sum_p[i][th2][0]/t1sum_p[th2]*(t1sum_p_err[th2]/t1sum_p[th2]+sum_p_err[i][th2][0]/sum_p[i][th2][0]));
 grate_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_pi[i][th2][0]/t1sum_pi[th2]);
 grate_pi_ac1[th1][th2]->SetPointError(i,emp[i],sum_pi[i][th2][0]/t1sum_pi[th2]*(t1sum_pi_err[th2]/t1sum_pi[th2]+sum_pi_err[i][th2][0]/sum_pi[i][th2][0]));
 //---- AC2 ----//
 grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]/t2sum_k[th1]);
 grate_k_ac2[th1][th2]->SetPointError(i,emp[i],sum_k[i][th1][1]/t2sum_k[th1]*(t2sum_k_err[th1]/t2sum_k[th1]+sum_k_err[i][th1][1]/sum_k[i][th1][1]));
 grate_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_p[i][th1][1]/t2sum_p[th1]);
 grate_p_ac2[th1][th2]->SetPointError(i,emp[i],sum_p[i][th1][1]/t2sum_p[th1]*(t2sum_p_err[th1]/t2sum_p[th1]+sum_p_err[i][th1][1]/sum_p[i][th1][1]));
 grate_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_pi[i][th1][1]/t2sum_pi[th1]);
 grate_pi_ac2[th1][th2]->SetPointError(i,emp[i],sum_pi[i][th1][1]/t2sum_pi[th1]*(t2sum_pi_err[th1]/t2sum_pi[th1]+sum_pi_err[i][th1][1]/sum_pi[i][th1][1]));
 */


 //----------- Rate TGraph ---------------------//
 //---- AC1 ----//
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]/sum_k[0][th2][0]);
 grate_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_p[i][th2][0]/sum_p[0][th2][0]);
 grate_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_pi[i][th2][0]/sum_pi[0][th2][0]);

 //---- AC2 ----//
 grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]/sum_k[0][th1][1]);
  grate_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_p[i][th1][1]/sum_p[0][th1][1]);
 grate_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_pi[i][th1][1]/sum_pi[0][th1][1]);






 //--------PointError -----------//

 gsum_k_ac1[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th2][0]);
 gsum_pi_ac1[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][0]);
 gsum_p_ac1[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th2][0]);
 gsum_pi_ac2[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][1]);
 gsum_k_ac2[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th1][1]);
 gsum_p_ac2[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th1][1]);
 grate_p_ac1[th1][th2]->SetPointError(i,0,sum_p[i][th2][0]/sum_p[0][th2][0]*(sum_p_err[i][th2][0]/sum_p[i][th2][0]+sum_p_err[0][th2][0]/sum_p[0][th2][0]));
 grate_k_ac1[th1][th2]->SetPointError(i,0,sum_k[i][th2][0]/sum_k[0][th2][0]*(sum_k_err[i][th2][0]/sum_k[i][th2][0]+sum_k_err[0][th2][0]/sum_k[0][th2][0]));
 grate_pi_ac1[th1][th2]->SetPointError(i,0,sum_pi[i][th2][0]/sum_pi[0][th2][0]*(sum_pi_err[i][th2][0]/sum_pi[i][th2][0]+sum_pi_err[0][th2][0]/sum_pi[0][th2][0]));
 grate_k_ac2[th1][th2]->SetPointError(i,emp[i],sum_k[i][th1][1]/sum_k[0][th1][1]*(sum_k_err[0][th1][1]/sum_k[0][th1][1]+sum_k_err[i][th1][1]/sum_k[i][th1][1]));
 grate_p_ac2[th1][th2]->SetPointError(i,emp[i],sum_p[i][th1][1]/t2sum_p[th1]*(sum_p_err[0][th1][1]/sum_p[0][th1][1]+sum_p_err[i][th1][1]/sum_p[i][th1][1]));
 grate_pi_ac2[th1][th2]->SetPointError(i,emp[i],sum_pi[i][th1][1]/t2sum_pi[th1]*(sum_pi_err[0][th1][1]/sum_pi[0][th1][1]+sum_pi_err[i][th1][1]/sum_pi[i][th1][1]));
//-----------------------------//

 if(n_k[i][th2][0]>max_nk[th1][th2][0])max_nk[th1][th2][0]=n_k[i][th2][0];   if(n_k[i][th1][1]>max_nk[th1][th2][1])max_nk[th1][th2][1]=n_k[i][th1][1];   
 if(n_p[i][th2][0]>max_np[th1][th2][0])max_np[th1][th2][0]=n_p[i][th2][0];   if(n_p[i][th1][1]>max_np[th1][th2][1])max_np[th1][th2][1]=n_p[i][th1][1];
 if(n_pi[i][th2][0]>max_npi[th1][th2][0])max_npi[th1][th2][0]=n_pi[i][th2][0];  if(n_pi[i][th1][1]>max_npi[th1][th2][1])max_npi[th1][th2][1]=n_pi[i][th1][1];
          
}

	}
      }

      cout<<" Fitting is done !"<<endl;
 
     


//======== Draw TCanvas ==============//




      TCanvas* ck1[nth];
      	for(int th2=0;th2<nth;th2++){
	  //	  if(th2==1){
 ck1[th2]=new TCanvas(Form("ck1[%d]",th2),Form("ck1[%d]",th2)); 
 ck1[th2]->cd();
     fk[0][th2][0]->SetLineColor(2);
     fk[0][th2][0]->SetFillColor(2);
     fk[0][th2][0]->SetFillStyle(3001);
     fpi[0][th2][0]->SetFillStyle(3001);
     fpi[0][th2][0]->SetLineColor(4);
     fpi[0][th2][0]->SetFillColor(4);
     fp[0][th2][0]->SetFillStyle(3001);
     fp[0][th2][0]->SetLineColor(3);
     fp[0][th2][0]->SetFillColor(3);
     hcoin_ac1_p[0][th2]->Draw();
     fk[0][th2][0]->Draw("same");
     fp[0][th2][0]->Draw("same");
     fpi[0][th2][0]->Draw("same");
     fcoin[0][th2][0]->Draw("same");
     //	}
      }



     
      TCanvas* ck2[nth];
       for(int th1=0;th1<nth;th1++){
	 //	  if(th1==1){
    ck2[th1]=new TCanvas(Form("ck2[%d]",th1),Form("ck2[%d]",th1)); 
     ck2[th1]->cd();
     fk[0][th1][1]->SetLineColor(2);
     fk[0][th1][1]->SetFillColor(2);
     fk[0][th1][1]->SetFillStyle(3001);
     fpi[0][th1][1]->SetFillStyle(3001);
     fpi[0][th1][1]->SetLineColor(4);
     fpi[0][th1][1]->SetFillColor(4);
     fp[0][th1][1]->SetFillStyle(3001);
     fp[0][th1][1]->SetLineColor(3);
     fp[0][th1][1]->SetFillColor(3);
     hcoin_ac2_p[0][th1]->Draw();
     fk[0][th1][1]->Draw("same");
     fp[0][th1][1]->Draw("same");
     fpi[0][th1][1]->Draw("same");
     fcoin[0][th1][1]->Draw("same");
     //	}
      }

     

 TLine* lac=new TLine(min_coin_c,ac1_adc[0],max_coin_c,ac1_adc[0]);
 lac->SetLineWidth(2);
 lac->SetLineColor(2);

 TCanvas* ccoin_ac=new TCanvas("ccoin_ac","ccoin_ac");
 ccoin_ac->Divide(1,2);
 ccoin_ac->cd(1); 
 hcoin_ac1[0]->Draw("colz");
 ccoin_ac->cd(2);
 hcoin_ac2[0]->Draw("colz");
 
 for(int j=0;j<nth;j++){
   ccoin_ac->cd(1);
 lac->SetLineColor(j+1);
 lac->DrawLine(min_coin_c,ac1_adc[j],max_coin_c,ac1_adc[j]);
 ccoin_ac->cd(2);
 lac->DrawLine(min_coin_c,ac2_adc[j],max_coin_c,ac2_adc[j]); 
 
}

 TCanvas* cAC=new TCanvas("cAC","cAC");
 cAC->cd();
 ha1_a2->Draw("colz");
 for(int j=0;j<nth;j++){
 lac->SetLineColor(j+1);
 lac->DrawLine(min_ac1,ac2_adc[j],max_ac1,ac2_adc[j]);
 lac->DrawLine(ac1_adc[j],min_ac2,ac1_adc[j],max_ac2);
 }



 lac->SetLineColor(2);
 
 /*
  TCanvas* ccoin_t1=new TCanvas("ccoin_t1","ccoin_t1");
  ccoin_t1->Divide(3,1);    
  int mv1=1;
  for(int th1=0;th1<nth;th1++){
    for(int th2=0;th2<nth;th2++){
      if(th1==1){
 fk_t1def[th1][th2]->SetLineColor(2);
 fp_t1def[th1][th2]->SetLineColor(3);
 fpi_t1def[th1][th2]->SetLineColor(4);
 fk_t1def[th1][th2]->SetFillStyle(3001);
 fp_t1def[th1][th2]->SetFillStyle(3001);
 fpi_t1def[th1][th2]->SetFillStyle(3001);
 fk_t1def[th1][th2]->SetFillColor(2);
 fp_t1def[th1][th2]->SetFillColor(3);
 fpi_t1def[th1][th2]->SetFillColor(4);
 
 ccoin_t1->cd(mv1);
      hcoin_t1[th2]->Draw();
      fcoin_t1[th1][th2]->Draw("same");
      fp_t1def[th1][th2]->Draw("same");   
      fpi_t1def[th1][th2]->Draw("same");
      fk_t1def[th1][th2]->Draw("same");
      mv1=mv1+1;
 }
  }
  }

  
  TCanvas* ccoin_t2=new TCanvas("ccoin_t2","ccoin_t2");
  ccoin_t2->Divide(3,1);    
  int mv2=1;
  for(int th1=0;th1<nth;th1++){
    for(int th2=0;th2<nth;th2++){
      if(th2==1){
 fk_t2def[th1][th2]->SetLineColor(2);
 fp_t2def[th1][th2]->SetLineColor(3);
 fpi_t2def[th1][th2]->SetLineColor(4);
 fk_t2def[th1][th2]->SetFillStyle(3001);
 fp_t2def[th1][th2]->SetFillStyle(3001);
 fpi_t2def[th1][th2]->SetFillStyle(3001);
 fk_t2def[th1][th2]->SetFillColor(2);
 fp_t2def[th1][th2]->SetFillColor(3);
 fpi_t2def[th1][th2]->SetFillColor(4);
 
 ccoin_t2->cd(mv2);
      hcoin_t2[th1]->Draw();
      fcoin_t2[th1][th2]->Draw("same");
      fp_t2def[th1][th2]->Draw("same");   
      fpi_t2def[th1][th2]->Draw("same");
      fk_t2def[th1][th2]->Draw("same");
      mv2=mv2+1;
 }
  }
  }

 */
 /*
     TCanvas* cAC1_cut=new TCanvas("cAC1_cut","cAC1_cut"); 
     cAC1_cut->Divide(4,3);
     for(int i=0;i<12;i++){
     for(int j=0;j<nth;j++){
     for(int k=0;k<nth;k++){
   int l=4*i;
       if(k==0 && j==0){
	 cAC1_cut->cd(i+1);
   
     fk[l][j][k][0]->SetLineColor(4);
     fk[l][j][k][0]->SetFillColor(4);
     fk[l][j][k][0]->SetFillStyle(3001);
     fpi[l][j][k][0]->SetFillStyle(3001);
     fpi[l][j][k][0]->SetLineColor(46);
     fpi[l][j][k][0]->SetFillColor(46);
     fp[l][j][k][0]->SetFillStyle(3001);
     fp[l][j][k][0]->SetLineColor(8);
     fp[l][j][k][0]->SetFillColor(8);
     hcoin_ac1_p[l][j][k]->Draw();
     fk[l][j][k][0]->Draw("same");
     fp[l][j][k][0]->Draw("same");
     fpi[l][j][k][0]->Draw("same");
     fcoin[l][j][k][0]->Draw("same");
       }
     }
    }
     }
 
 

  TCanvas* cAC2_cut=new TCanvas("cAC2_cut","cAC2_cut"); 
     cAC2_cut->Divide(4,3);
  for(int i=0;i<12;i++){
     for(int j=0;j<nth;j++){
     for(int k=0;k<nth;k++){
       if(j==0 && k==0){
   int  l=4*i;
   cAC2_cut->cd(i+1);
     hcoin_ac2_p[l][j][k]->Draw();
     fk[l][j][k][1]->SetLineColor(4);
     fk[l][j][k][1]->SetFillColor(4);
     fk[l][j][k][1]->SetFillStyle(3001);
     fpi[l][j][k][1]->SetFillStyle(3001);
     fpi[l][j][k][1]->SetLineColor(46);
     fpi[l][j][k][1]->SetFillColor(46);
     fp[l][j][k][1]->SetFillStyle(3001);
     fp[l][j][k][1]->SetLineColor(8);
     fp[l][j][k][1]->SetFillColor(8);
     fk[l][j][k][1]->Draw("same");
     fp[l][j][k][1]->Draw("same");
     fpi[l][j][k][1]->Draw("same");
     fcoin[l][j][k][1]->Draw("same");
     

       }
     }
     }
  }

   





     
 */

    
  
 TCanvas* cAC1_num=new TCanvas("cAC1_num","cAC1_num");
          cAC1_num->Divide(1,3);

 TCanvas* cAC1_rate=new TCanvas("cAC1_rate","cAC1_rate");
          cAC1_rate->Divide(1,3);
      	  
       for(int k=0;k<nth;k++){
    for(int j=0;j<nth;j++){
      if(j==1){
       if(k==0){
   cAC1_num->cd(1);
   gsum_k_ac1[j][k]->SetMarkerColor(1);
   gsum_k_ac1[j][k]->SetFillColor(1);
   gsum_k_ac1[j][k]->SetFillStyle(3005);
   gsum_k_ac1[j][k]->Draw("Ap");
   cAC1_num->cd(2);
   gsum_pi_ac1[j][k]->SetFillColor(1);
   gsum_pi_ac1[j][k]->SetFillStyle(3005);
   gsum_pi_ac1[j][k]->SetMarkerColor(1);
   gsum_pi_ac1[j][k]->Draw("Ap");
   cAC1_num->cd(3);
   gsum_p_ac1[j][k]->SetFillColor(1);
   gsum_p_ac1[j][k]->SetFillStyle(3005);
   gsum_p_ac1[j][k]->SetMarkerColor(1);
   gsum_p_ac1[j][k]->Draw("Ap");
   
    cAC1_rate->cd(1);
    // grate_k_ac1[j][k]->SetFillColor(1);
   // grate_k_ac1[j][k]->SetFillStyle(3005);
   grate_k_ac1[j][k]->SetMarkerColor(1);
   grate_k_ac1[j][k]->Draw("Ap");
   cAC1_rate->cd(2);
   // grate_pi_ac1[j][k]->SetFillColor(1);
   // grate_pi_ac1[j][k]->SetFillStyle(3005);
   grate_pi_ac1[j][k]->SetMarkerColor(1); 
   grate_pi_ac1[j][k]->Draw("Ap");
   cAC1_rate->cd(3);
   // grate_k_ac1[j][k]->SetFillColor(1);
   // grate_k_ac1[j][k]->SetFillStyle(3005);
   grate_p_ac1[j][k]->SetMarkerColor(1);
   grate_p_ac1[j][k]->Draw("Ap");
     
   
      }else{
   cAC1_num->cd(1);
   gsum_k_ac1[j][k]->SetMarkerColor(k+1);
   gsum_k_ac1[j][k]->Draw("P");
   cAC1_num->cd(2);   
   gsum_pi_ac1[j][k]->SetMarkerColor(k+1);
   gsum_pi_ac1[j][k]->Draw("P");
   cAC1_num->cd(3);   
   gsum_p_ac1[j][k]->SetMarkerColor(k+1);
   gsum_p_ac1[j][k]->Draw("P");
   
   cAC1_rate->cd(1);
   grate_k_ac1[j][k]->SetMarkerColor(k+1);
   grate_k_ac1[j][k]->Draw("P");
   cAC1_rate->cd(2);   
   grate_pi_ac1[j][k]->SetMarkerColor(k+1);
   grate_pi_ac1[j][k]->Draw("P");
   cAC1_rate->cd(3);   
   grate_p_ac1[j][k]->SetMarkerColor(k+1);
   grate_p_ac1[j][k]->Draw("P");
    }
       }
    }
       }

 TCanvas* cAC2_rate=new TCanvas("cAC2_rate","cAC2_rate");
          cAC2_rate->Divide(1,3);

 TCanvas* cAC2_num=new TCanvas("cAC2_num","cAC2_num");
          cAC2_num->Divide(1,3);
	         for(int j=0;j<nth;j++){
       for(int k=0;k<nth;k++){
	 if(k==1){
	  if(j==0){
   cAC2_num->cd(1);
   gsum_k_ac2[j][k]->SetMarkerColor(1);
   gsum_k_ac2[j][k]->Draw("AP");
   cAC2_num->cd(2);
   gsum_pi_ac2[j][k]->SetMarkerColor(1);
   gsum_pi_ac2[j][k]->Draw("AP");
   cAC2_num->cd(3);
   gsum_p_ac2[j][k]->SetMarkerColor(1);
   gsum_p_ac2[j][k]->Draw("AP");

   cAC2_rate->cd(1);
   grate_k_ac2[j][k]->SetMarkerColor(1);
   grate_k_ac2[j][k]->Draw("AP");
   cAC2_rate->cd(2);
   grate_pi_ac2[j][k]->SetMarkerColor(1);
   grate_pi_ac2[j][k]->Draw("AP");
   cAC2_rate->cd(3);
   grate_p_ac2[j][k]->SetMarkerColor(1);
   grate_p_ac2[j][k]->Draw("AP");
	  }else{
   cAC2_num->cd(1);
   gsum_k_ac2[j][k]->SetMarkerColor(j+1);
   gsum_k_ac2[j][k]->Draw("P");
   cAC2_num->cd(2);
   gsum_pi_ac2[j][k]->SetMarkerColor(j+1);
   gsum_pi_ac2[j][k]->Draw("P");
   cAC2_num->cd(3);
   gsum_p_ac2[j][k]->SetMarkerColor(j+1);
   gsum_p_ac2[j][k]->Draw("P");

   cAC2_rate->cd(1);
   grate_k_ac2[j][k]->SetMarkerColor(j+1);
   grate_k_ac2[j][k]->Draw("P");
   cAC2_rate->cd(2);
   grate_pi_ac2[j][k]->SetMarkerColor(j+1);
   grate_pi_ac2[j][k]->Draw("P");
   cAC2_rate->cd(3);
   grate_p_ac2[j][k]->SetMarkerColor(j+1);
   grate_p_ac2[j][k]->Draw("P");

	  }
	 }
       }
        }

 
      cout<<"drawing is done !!"<<endl;       

 //================ Print Canvas =================================//

        
TString name;
 if(output_flag){
 
 name.Form(ofname.c_str());
 ccoin_ac->Print(name+"[","pdf");
 ccoin_ac->Print(name,"pdf");
 cAC->Print(name,"pdf"); 


 for(int j=0;j<nth;j++){
ck1[j]->Print(name,"pdf");
ck2[j]->Print(name,"pdf");
 }

 cAC1_num->Print(name,"pdf");
 cAC2_num->Print(name,"pdf");
 cAC1_rate->Print(name,"pdf");
 cAC2_rate->Print(name,"pdf");
 cAC2_rate->Print(name+"]","pdf");

 

 // cAC1_cut->Print(name,"pdf");
 //cAC2_cut->Print(name,"pdf");
 //cAC2_cut->Print(name +"]","pdf");
 // cAC1_cut2->Print(name,"pdf");
 // cAC1_cut3->Print(name,"pdf");
 // cAC1_cut4->Print(name,"pdf");
 // cAC2_cut2->Print(name,"pdf");
 // cAC2_cut3->Print(name,"pdf");
 // cAC2_cut4->Print(name,"pdf");
 // cAC2_cut4->Print(name +"]","pdf");


 }     
    
 cout<<"Print is done "<<endl;
   


 for(int j=0;j<nth;j++){
   for(int k=0;k<nth;k++){

 cout<<"========= j:"<<j<<" k:"<<k<<"========================"<<endl; 

 cout<<Form("ac1_adc[%d]:",j)<<ac1_adc[j]<<Form(" : ac2_adc[%d]:",k)<<ac2_adc[k]<<endl;
 cout<<"max k[ac1]: "<<max_nk[j][k][0]/tdc_time<<" : max k[ac2]: "<<max_nk[j][k][1]/tdc_time<<endl;

 cout<<"========================================================="<<endl;
   }
 }
 
 if(draw_flag==0)gSystem->Exit(1);
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
