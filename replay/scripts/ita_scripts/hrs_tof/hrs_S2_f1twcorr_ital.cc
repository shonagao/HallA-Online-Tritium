// 10/4/2018 Author Itabshi
// time walk correction R-HRS S0 F1TDC

// Grobal Function //
 int chmax=17; // channel of S2 PMT 
// int ch=8;// defolt ch 8;
 double tdc_time=56.25e-12;//TDC converse ch->sec [sec/ch]
int ch=8;
bool s2=true;

//====== Define Function ===========//
double range_para(int i,int j);
double range_adc(int i,int j);
long double fit_ini(int i,int j,int k);
//=================================//
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

//-------- TTree data input ---------------//
  TChain*  T=new TChain("T");


  T->Add(Form("/w/halla-scifs17exp/triton/itabashi/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",94003));

 gStyle->SetOptFit(0001);
 //=================================================//
 //========== Get Tree Branch ======================//
 //=================================================//

  int max=10000; 
  double F1[max];
  double s0radc;
  double s0ladc;
  double s2radc[max];
  double s2ladc[max];
  int trig;
  bool rarm=true;
 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",F1); 
 
  if(rarm){
  T->SetBranchStatus("R.s0.ra_c",1);        // Right arm S0 R-PMT  ADC
  T->SetBranchAddress("R.s0.ra_c",&s0radc); // Right arm S0 R-PMT  ADC
  T->SetBranchStatus("R.s0.la_c",1);        // Right arm S0 L-PMT  ADC
  T->SetBranchAddress("R.s0.la_c",&s0ladc); // Right arm S0 L-PMT  ADC
  T->SetBranchStatus("R.s2.ra_c",1);        // Right arm S2 R-PMT  ADC
  T->SetBranchAddress("R.s2.ra_c",s2radc);  // Right arm S2 R-PMT  ADC
  T->SetBranchStatus("R.s2.la_c",1);        // Right arm S2 L-PMT  ADC
  T->SetBranchAddress("R.s2.la_c",s2ladc);  // Right arm S2 L-PMT  ADC
  }else{
  T->SetBranchStatus("L.s0.ra_c",1);        // Left arm S0 R-PMT  ADC
  T->SetBranchAddress("L.s0.ra_c",&s0radc); // Left arm S0 R-PMT  ADC
  T->SetBranchStatus("L.s0.la_c",1);        // Left arm S0 L-PMT  ADC
  T->SetBranchAddress("L.s0.la_c",&s0ladc); // Left arm S0 L-PMT  ADC
  T->SetBranchStatus("L.s2.ra_c",1);        // Left arm S2 R-PMT  ADC
  T->SetBranchAddress("L.s2.ra_c",s2radc);  // Left arm S2 R-PMT  ADC
  T->SetBranchStatus("L.s2.la_c",1);        // Left arm S2 L-PMT  ADC
  T->SetBranchAddress("L.s2.la_c",s2ladc);  // Left arm S2 L-PMT  ADC
  }
 

 //=================================================//
 //=============== Definition ======================//
 //=================================================//

  //==== w/o Time-Walk Correction Paramters =====//
  //--- TCanvas ----//
  TCanvas* ctof[chmax];
  TCanvas* cwo[chmax];
  //--- Function ----//
  TH1F* htof[chmax];
  TH2F* htof_adc_r[chmax];
  TH2F* htof_adc_l[chmax];
  //--- Parameters ---//
  int evnt;
  double tof[chmax];
  double bin_tof[chmax],min_tof[chmax],max_tof[chmax];
  double bin_adc[chmax],min_adc[chmax],max_adc[chmax];
  double bin_tof_c[chmax],min_tof_c[chmax],max_tof_c[chmax];
  double bin_adc_c[chmax],min_adc_c[chmax],max_adc_c[chmax];
  double S0a_t, S0b_t,S2l_t[chmax],S2r_t[chmax];
  double L1,L2;
  double S0a_a,S0b_a,S2l_a[chmax],S2r_a[chmax];
  // Cut conditon //
  bool tdc_cut[chmax];
  bool adc_cut[chmax];
  bool trig_cut;
  bool tof_cut[chmax];
  //======= TSlicesY Paramters =======//
  TCanvas* cslice[chmax];
  TH1F* hslice_tof_r[chmax];
  TH1F* hslice_tof_l[chmax];
  TF1* ftof_gaus[chmax];
  TF1* ffit_r[chmax];
  TF1* fcal[chmax];
  TF1* ffit_l[chmax];
  int tw_nparam=2;
  double twp_r[chmax][tw_nparam],twp_l[chmax][tw_nparam];
 //====== Fill w/ Time-Walk Correction ======//
  TCanvas* cw[chmax];
  TH1F* htof_c[chmax];
  TH2F *htof_adc_rc[chmax]; 
  TH2F *htof_adc_lc[chmax];
  double tof_c[chmax];
  double corr_r[chmax],corr_l[chmax];


  //=========== Definition Functions =============================//

  for(int i=0;i<chmax;i++){
    if(i==ch){
      //------------ w/o Time-Walk Correction Parameters ----------------------//
      // min_tof[i]=-2.5e-8;   max_tof[i]=-1.5e-8;//1.0e-6;
       min_tof[i]=0.0;   max_tof[i]=3.0e-7;//1.0e-6;
      bin_tof[i]=(max_tof[i]-min_tof[i])/tdc_time; bin_tof[i]=(int)bin_tof[i];
      min_adc[i]=0.0;  max_adc[i]=500.;
      bin_adc[i]=(max_adc[i]-min_adc[i]);      bin_adc[i]=(int)bin_adc[i];
      ctof[i]=new TCanvas(Form("ctof[%d]",i),Form("ctof[%d]",i));
      ctof[i]->Divide(1,2);
      cwo[i]=new TCanvas(Form("cwo[%d]",i),Form("cwo[%d]",i)); 
      cwo[i]->Divide(1,2);
      htof[i]=new TH1F(Form("htof[%d]",i),Form("S2ch%d -S0 TOF w/o Time-Walk Correction",i),bin_tof[i],min_tof[i],max_tof[i]);
      htof_adc_r[i]=new TH2F(Form("htof_adc_r[%d]",i),Form("S2ch%d -S0 TOF vs S2 ch%d R-PMT ADC w/o Time-Walk Correction",i,i),bin_adc[i],min_adc[i],max_adc[i],bin_tof[i],min_tof[i],max_tof[i]);
      htof_adc_l[i]=new TH2F(Form("htof_adc_l[%d]",i),Form("S2ch%d -S0 TOF vs S2 ch%d L-PMT ADC w/o Time-Walk Correction",i,i),bin_adc[i],min_adc[i],max_adc[i],bin_tof[i],min_tof[i],max_tof[i]);

      //------------ TSlicesY Parameters ----------------------//
      cslice[i]=new TCanvas(Form("cslice[%d]",i),Form("cslice[%d]",i));
      ftof_gaus[i]=new TF1(Form("ftof_gaus[%d]",i),"gaus",min_tof[i],max_tof[i]);
      ffit_r[i]=new TF1(Form("ffit_r[%d]",i),"[0]*1./sqrt(x)-[1]",min_adc[i],max_adc[i]);
      ffit_l[i]=new TF1(Form("ffit_l[%d]",i),"[0]*1./sqrt(x)-[1]",min_adc[i],max_adc[i]);      

      //------------ w/ Time-Walk Correction Parameters ----------------------//
      min_tof_c[i]=1.0e-8;   max_tof_c[i]=6.0e-8;//1.0e-6;
      bin_tof_c[i]=(max_tof_c[i]-min_tof_c[i])/tdc_time;  bin_tof_c[i]=(int)bin_tof_c[i];
      min_adc_c[i]=100.0;   max_adc_c[i]=1000.;
      bin_adc_c[i]=max_adc_c[i]-min_adc_c[i];  bin_adc_c[i]=(int)bin_adc_c[i];
    htof_c[i]=new TH1F(Form("htof_c[%d]",i),Form("S2ch%d -S0 TOF w/ Time-Walk Correction",i),bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
    htof_adc_rc[i]=new TH2F(Form("htof_adc_rc[%d]",i),Form("S2ch%d -S0 TOF vs S2 ch%d R-PMT ADC w/ Time-Walk Correction",i,i),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
      htof_adc_lc[i]=new TH2F(Form("htof_adc_lc[%d]",i),Form("S2ch%d -S0 TOF vs S2 ch%d L-PMT ADC w/ Time-Walk Correction",i,i),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);

 cw[i]=new TCanvas(Form("cw[%d]",i),Form("cw[%d]",i)); 
     
    }
  }

 //=========================================================//
 


 //===============================================================//
 //========== Fill w/o Time-Walk Correciton ======================//
 //===============================================================//

  evnt=T->GetEntries();
  cout<<"Get Entry : "<<evnt<<endl;
  for(int i=0;i<chmax;i++){
    if(i==ch){
      for(int k=0;k<evnt;k++){
	T->GetEntry(k);
 
 //---- Cut condition --//
  tdc_cut[i]=false;
  adc_cut[i]=false;
  trig_cut=false;
  tof_cut[i]=false;
 if(rarm){
  S0a_t=(-F1[43]+F1[46])*tdc_time;
  S0b_t=(-F1[44]+F1[46])*tdc_time;
  S2l_t[i]=(-F1[16+i]+F1[9])*tdc_time;
  S2r_t[i]=(-F1[48+i]+F1[46])*tdc_time;
  if(F1[43]>0 && F1[44]>0 && F1[16+i]>0 && F1[48+i]>0)tdc_cut[i]=true;
  if(F1[4]>0)trig_cut=true; 
  if(tof[i]>1.6e-8 && tof[i]<3.0e-8)tof_cut[i]=true;
 }else{
  S0a_t=(-F1[27]+F1[30])*tdc_time;
  S0b_t=(-F1[28]+F1[30])*tdc_time;
  S2l_t[i]=(-F1[0+i]+F1[30])*tdc_time;
  S2r_t[i]=(-F1[48+i]+F1[40])*tdc_time;
  if(F1[27]>0 && F1[28]>0 && F1[0+i]>0 && F1[48+i]>0)tdc_cut[i]=true;
  if(F1[33]>0)trig_cut=true; 
  if(tof[i]>1.6e-8 && tof[i]<2.6e-8)tof_cut[i]=true;
  }
 
  S0a_a=s0radc;
  S0b_a=s0ladc;
  S2r_a[i]=s2radc[i];
  S2l_a[i]=s2ladc[i];
  tof[i]=(S2l_t[i]+S2r_t[i])/2.0-(S0a_t+S0b_t)/2.0;
 
  //--------- Fill Hist ---------------//  
  if(tdc_cut[i] && trig_cut && tof_cut[i])htof[i]->Fill(tof[i]);
  if(tdc_cut[i] && trig_cut&& tof_cut[i])htof_adc_r[i]->Fill(S2r_a[i],tof[i]);
  if(tdc_cut[i] && trig_cut&& tof_cut[i])htof_adc_l[i]->Fill(S2l_a[i],tof[i]);

      }
    
    }
  }




 //=====================================================================//
 //========== TSlicesY Get Correction Parameters ======================//
 //====================================================================//



    for(int i=0;i<chmax;i++){  
     if(i==ch){
 //------ FitSlicesY -----------------// 
 TObjArray aSlices;     
 htof_adc_r[i]->FitSlicesY(ftof_gaus[i],0,-1,0,"QRG3");
 hslice_tof_r[i]=(TH1F*)gROOT->FindObject(Form("htof_adc_r[%d]_1",i));
 htof_adc_l[i]->FitSlicesY(ftof_gaus[i],0,-1,0,"QRG3");
 hslice_tof_l[i]=(TH1F*)gROOT->FindObject(Form("htof_adc_l[%d]_1",i));

     
 //===== S0 Set Paramters =======//
 //#94003 run//
 if(rarm){
  ffit_r[i]->SetParameter(0,7.3e-8);   
  ffit_r[i]->SetParameter(1,-1.87e-8);   
  ffit_l[i]->SetParameter(0,1.0e-7);   
  ffit_l[i]->SetParameter(1,1.0-8);    
  ffit_l[i]->SetParLimits(0,1.0e-7,4.0e-7);   
 }else{
  ffit_r[i]->SetParameter(0,7.3e-8);   
  ffit_r[i]->SetParameter(1,-1.87e-8);   
  ffit_l[i]->SetParameter(0,1.0e-7);   
  ffit_l[i]->SetParameter(1,1.0-8);    
  ffit_l[i]->SetParLimits(0,1.0e-7,4.0e-7);   
}
  //==============================//
 hslice_tof_r[i]->Fit(Form("ffit_r[%d]",i),"","");
 hslice_tof_l[i]->Fit(Form("ffit_l[%d]",i),"","");


      twp_r[i][0]=ffit_r[i]->GetParameter(0);
      twp_r[i][1]=ffit_r[i]->GetParameter(1);
      twp_l[i][0]=ffit_l[i]->GetParameter(0);
      twp_l[i][1]=ffit_l[i]->GetParameter(1);
     }
    }  
    


 //===============================================================//
 //=========== Fill w/ Time-Walk Correciton ======================//
 //===============================================================//

    //itallation parameters //
    int ital_a=1;
    int ital_b=1;
    int amin[chmax],bmin[chmax];
    double wa[chmax],wb[chmax],wamin[chmax],wbmin[chmax];
    double sig_ital[chmax][ital_a][ital_b],sig_c[chmax];
    for(int i=0;i<chmax;i++){ sig_c[i]=100.;}
    TH1F* htof_ital[chmax][ital_a][ital_b];
    TF1* ftof_c[chmax];
    TH2F* htof_adc_italr[chmax][ital_a][ital_b];
    TH2F* htof_adc_itall[chmax][ital_a][ital_b];
    TH2F* hital[chmax];

    for(int i=0;i<chmax;i++){ 
     if(i==ch){
 //------------ TOF italtion  hist -------------//
       hital[i]=new TH2F(Form("hital[%d]",i),Form("TOF Itallation hist ch %d",i),ital_a,0,ital_a,ital_b,0,ital_b);
ftof_c[i]=new TF1(Form("ftof_c[%d]",i),"gaus",min_tof_c[i],max_tof_c[i]);
       for(int a=0;a<ital_a;a++){
        for(int b=0;b<ital_b;b++){

 htof_ital[i][a][b]=new TH1F(Form("htof_ital[%d][%d][%d]",i,a,b),Form("TOF with Scinti CH %d, w/ correlation ",i),bin_tof_c[i],min_tof_c[i],max_tof_c[i]); 
  //---- TOF vs ADC 2D Hist w/ Time-Walk Correction ---//
 htof_adc_italr[i][a][b]=new TH2F(Form("htof_adc_italr[%d][%d][%d]",i,a,b),Form("TOF vs Scinti CH %d,R-PMT adc w/ correlation ",i),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
htof_adc_itall[i][a][b]=new TH2F(Form("htof_adc_itall[%d][%d][%d]",i,a,b),Form("TOF vs Scinti CH %d,L-PMT adc w/ correlation ",i),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);


         for(int k=0;k<evnt;k++){ // Fill Entr
        T->GetEntry(k);

 //---- Cut condition --//
  tdc_cut[i]=false;
  adc_cut[i]=false;
  trig_cut=false;
  tof_cut[i]=true;
 //----- F1TDC ------// 
 if(rarm){
  S0a_t=(-F1[43]+F1[46])*tdc_time;
  S0b_t=(-F1[44]+F1[46])*tdc_time;
  S2l_t[i]=(-F1[16+i]+F1[9])*tdc_time;
  S2r_t[i]=(-F1[48+i]+F1[46])*tdc_time;
  if(F1[43]>0 && F1[44]>0 && F1[16+i]>0 && F1[48+i]>0)tdc_cut[i]=true;
  if(F1[4]>0)trig_cut=true; 
 if(tof[i]>1.6e-8 && tof[i]<2.6e-8)tof_cut[i]=true;
 }else{
  S0a_t=(-F1[27]+F1[30])*tdc_time;
  S0b_t=(-F1[28]+F1[30])*tdc_time;
  S2l_t[i]=(-F1[0+i]+F1[30])*tdc_time;
  S2r_t[i]=(-F1[48+i]+F1[40])*tdc_time;
  if(F1[27]>0 && F1[28]>0 && F1[0+i]>0 && F1[48+i]>0)tdc_cut[i]=true;
  if(F1[33]>0)trig_cut=true; 
 if(tof[i]>1.6e-8 && tof[i]<2.6e-8)tof_cut[i]=true;
  }

 //------FBUS ADC ----//
  S0a_a=s0radc;
  S0b_a=s0ladc;
  S2r_a[i]=s2radc[i];
  S2l_a[i]=s2ladc[i];
  //---- Time-Walk Correction Parameters ---//
  //  corr_r[i]=twp_r[i][0]*1./S0a_a-twp_r[i][1];
  // corr_l[i]=twp_l[i][0]*1./S0b_a-twp_l[i][1];
  corr_r[i]=twp_r[i][0]*1./S0a_a;
  corr_l[i]=twp_l[i][0]*1./S0b_a;
	wa[i]=0.0+0.1*a;
	wb[i]=0.0+0.1*b;
	tof[i]=(S2l_t[i]+S2r_t[i])/2.0-(S0a_t+S0b_t)/2.0;
  tof_c[i]=(S2l_t[i]+S2r_t[i])/2.0-(S0a_t+S0b_t)/2.0-wa[i]*corr_r[i]-wb[i]*corr_l[i];
  //------- Fill Hist w/ Time-Walk Correction ------//

  if(tdc_cut[i] && trig_cut && tof_cut[i])htof_ital[i][a][b]->Fill(tof_c[i]);
  if(tdc_cut[i] && trig_cut && tof_cut[i])htof_adc_italr[i][a][b]->Fill(S2r_a[i],tof_c[i]);
  if(tdc_cut[i] && trig_cut && tof_cut[i])htof_adc_itall[i][a][b]->Fill(S2l_a[i],tof_c[i]);

      } //END Fill


 htof_ital[i][a][b]->Fit(ftof_c[i]);
 sig_ital[i][a][b]=ftof_c[i]->GetParameter(2);
 if(sig_ital[i][a][b]>1.0e-10){
 if(sig_c[i]>sig_ital[i][a][b]){
 sig_c[i]=sig_ital[i][a][b]; 
 amin[i]=a; bmin[i]=b;
 wamin[i]=wa[i];wbmin[i]=wb[i];
 htof_c[i]=(TH1F*)htof_ital[i][a][b]->Clone(); 
}
 }
 cout<<"i "<<i<<": ital_a "<<a+1<<Form("/%d :ital_b ",ital_a)<<b+1<<Form("/%d",ital_b)<<endl;
 cout<<"sig_ital :"<<sig_ital[i][a][b]<<endl;
 cout<<"sig_c :"<<sig_c[i]<<endl;
 cout<<"amin "<<amin[i]<<" : bmin "<<bmin[i]<<endl;
 cout<<"wamin "<<wamin[i]<<" : wbmin "<<wbmin[i]<<endl;
 hital[i]->Fill(a,b,sig_ital[i][a][b]);

	}
       }


     }
    }
    

    //=========================================================//
    //=============< Draw Hist >===============================// 
    //=========================================================//

    TCanvas* cital[chmax];
    bool draw=true;
    for(int i=0;i<chmax;i++){ 
     if(i==ch && draw){
    //---- w/o Time-Walk Correction Hist -----// 
      ctof[i]->cd(1);
      htof[i]->Draw();
      cwo[i]->cd(1);
      htof_adc_r[i]->Draw();
      cwo[i]->cd(2);
      htof_adc_l[i]->Draw();
    //------------ TSlicesY Hist ------------//
      cslice[i]->Divide(1,2);
      cslice[i]->cd(1);
      htof_adc_r[i]->Draw();
      hslice_tof_r[i]->Draw("same");
      ffit_r[i]->Draw("same");
      cslice[i]->cd(2);
      htof_adc_l[i]->Draw();
      hslice_tof_l[i]->Draw("same");
      ffit_l[i]->Draw("same");
   //---- w/ Time-Walk Correction Hist -----// 
      cw[i]->Divide(1,2);
      cw[i]->cd(1);
      // htof_adc_rc[i]->Draw();
      htof_adc_italr[i][amin[i]][bmin[i]]->Draw();
      cw[i]->cd(2);
       htof_adc_itall[i][amin[i]][bmin[i]]->Draw();
      //htof_adc_rc[i]->Draw();
      ctof[i]->cd(2);
      // htof_c[i]->Draw();
      htof_ital[i][amin[i]][bmin[i]]->Draw();
//========= TOF Itallation Hist =========//
   cital[i]=new TCanvas(Form("cital[%d]",i),Form("cital[%d]",i));
   cital[i]->cd(1);
   hital[i]->Draw("colz");

    }
    }
    
    cout<<"Draw is done "<<endl;

 theApp->Run();
 return 0;


}// END Main




//=======================================================================//
//=================== Definiton Function ================================//
//=======================================================================//



double range_para(int i,int j){

  int npara=6;// number of parameters
  double par[chmax][npara];
  double param;
 
  //=== Inital parameters========//
  for(int k=0;k<chmax;k++){
   
   par[k][0]=-4e-9, par[k][1]=4e-9, par[k][2]=-0.1e-7, par[k][3]=0.1e-7;//TOF

   if(k==16){//S0
  par[k][4]=500.0, par[k][5]=5000.;//ADC
    }else{ //S2
      par[k][4]=100.0, par[k][5]=500.;//ADC
    }
   }

  //===== Set Parameters ========//
  //  par[8][0]=0.18e-6, par[8][1]=0.22e-6, par[8][2]=-0.18e-6,par[8][3]=0.18e-6;
  //  par[16][0]=0.180e-6, par[16][1]=0.23e-6, par[16][2]=-0.18e-6, par[16][3]=0.25e-6;

  //============================//
  for(int k=0;k<chmax;k++){
   for(int l=0;l<npara;l++){
    if(k==i && l==j){
      param=par[k][l];
     
    }

    }
  }
  return param;
}






long double fit_ini(int i,int j,int k){
  // COMMENT  i: ch num, j: R-PMT(0),L-PMT(1), k: para num 
  int npara=2;
  double fpar[chmax][2][npara];
  double f_ini_param;

  //======== Initial parameters =========//

  for(int k=0;k<chmax;k++){
  fpar[k][0][0]= -6.0e-9, fpar[k][0][1]=1.0e-3,  fpar[k][1][0]= -1.25e-8, fpar[k][1][1]=3.67e-3;
   }



  // 2018 Oct. 20th
  fpar[7][0][0]= -1.38e-8, fpar[7][0][1]=1.0e5,  fpar[7][1][0]= -1.25e-7, fpar[7][1][1]=3.67e5;
  fpar[8][0][0]=-2.0e-8;


  fpar[8][0][0]= -5.7e-7, fpar[8][0][1]=0.81e5,  fpar[8][1][0]= -1.25e-7, fpar[8][1][1]=3.67e5;
  fpar[8][0][0]=-2.0e-8;

  //2018 Oct. 17th
 fpar[8][0][0]= -1.67e-8, fpar[8][0][1]=3.59e2,  fpar[8][1][0]= -1.85e-8, fpar[8][1][1]=3.34e2;
 

  //fpar[16][0][0]= -2.55e-7, fpar[16][0][1]=1.61,  fpar[16][1][0]= -7.40e-7, fpar[16][1][1]=1.0e2;
fpar[16][0][0]= -2.55e-9, fpar[16][0][1]=5.0e3,  fpar[16][1][0]= -3.40e-8, fpar[16][1][1]=5.0e3;
//2018 Oct. 17th
fpar[16][0][0]= -1.08e-8, fpar[16][0][1]=8.17e13,  fpar[16][1][0]= -1.11e-8, fpar[16][1][1]=1.92e15;

  for(int x=0;x<chmax;x++){
    for(int y=0;y<2;y++){
      for(int z=0;z<npara;z++){
       if(x==i && y==j && z==k){
        f_ini_param=fpar[i][j][k];
      }
     }
    }
  }
  return f_ini_param;
}


