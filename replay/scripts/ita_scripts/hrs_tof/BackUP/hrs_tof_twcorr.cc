// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2


// Grobal Function //
 int chmax=17; // channel of S2 PMT 
// int ch=8;// defolt ch 8;
 double tdc_time=0.1e-9;//TDC converse ch->sec [ss/ch]
// Itallation parameters
int ital1=10;
int ital2=10;
////////////////////////////
//Channel Select S2 ch 0 to 15//
int ch_st=0;
int ch_end=15;
///////////////////////////
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
  TChain*  T=new TChain("T");;
 for(Int_t i = 93495; i<93504;i++) /// these run goes 93495 to 93504
 {
    T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }
 gStyle->SetOptFit(0001);
 //=================================================//
 //========== Get Tree Branch ======================//
 //=================================================//

int max=10000; 
 int Ndata_Rs0la_c,Ndata_Rs0ra_c,Ndata_Rs0rt_c,Ndata_Rs0lt_c;//S0 Ndata
 int Ndata_Rs2la_c,Ndata_Rs2ra_c,Ndata_Rs2rt_c,Ndata_Rs2lt_c;//S2 Ndata
 double DR_evtypebits; // trigger condition  
 double Rs0la_c[max],Rs0ra_c[max],Rs0lt_c[max],Rs0rt_c[max]; //S0 TDC & ADC  
 double Rs2la_c[max],Rs2ra_c[max],Rs2lt_c[max],Rs2rt_c[max]; //S2 TDC & ADC 

 
 //============= Set Branch Status ==================//
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.R.s0.la_c",1);
  T->SetBranchStatus("Ndata.R.s0.ra_c",1);
  T->SetBranchStatus("Ndata.R.s0.lt_c",1);
  T->SetBranchStatus("Ndata.R.s0.rt_c",1);
  T->SetBranchStatus("Ndata.R.s2.la_c",1);
  T->SetBranchStatus("Ndata.R.s2.ra_c",1);
  T->SetBranchStatus("Ndata.R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1); 
  T->SetBranchStatus("R.s0.la_c",1);
  T->SetBranchStatus("R.s0.ra_c",1);
  T->SetBranchStatus("R.s0.lt_c",1);
  T->SetBranchStatus("R.s0.rt_c",1);
  T->SetBranchStatus("R.s2.la_c",1);
  T->SetBranchStatus("R.s2.ra_c",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("DR.evtypebits",1);   

  // Right Arm S0 Ndata// 
  T->SetBranchAddress("Ndata.R.s0.la_c",&Ndata_Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s0.ra_c",&Ndata_Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("Ndata.R.s2.la_c",&Ndata_Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("Ndata.R.s2.ra_c",&Ndata_Rs2ra_c); // Right arm S2-Bottom(A) ADC
  // Right Arm S2 Ndata//
  T->SetBranchAddress("Ndata.R.s0.lt_c",&Ndata_Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s0.rt_c",&Ndata_Rs0rt_c); // Right arm S0-Bottom(A) TDC
  T->SetBranchAddress("Ndata.R.s2.lt_c",&Ndata_Rs2lt_c); // Right arm S2-Top(B) TDC
  T->SetBranchAddress("Ndata.R.s2.rt_c",&Ndata_Rs2rt_c); // Right arm S2-Bottom(A) TDC
  // Right Arm S0 // 
  T->SetBranchAddress("R.s0.la_c",Rs0la_c); // Right arm S0-Top(B) ADC
  T->SetBranchAddress("R.s0.ra_c",Rs0ra_c); // Right arm S0-Bottom(A) ADC
  T->SetBranchAddress("R.s0.lt_c",Rs0lt_c); // Right arm S0-Top(B) TDC
  T->SetBranchAddress("R.s0.rt_c",Rs0rt_c); // Right arm S0-Bottom(A) TDC
  // Right Arm S2 // 
  T->SetBranchAddress("R.s2.la_c",Rs2la_c); // Right arm S2-Top(B) ADC
  T->SetBranchAddress("R.s2.ra_c",Rs2ra_c); // Right arm S2-Bottom(A) ADC
  T->SetBranchAddress("R.s2.lt_c",Rs2lt_c); // Right arm Se2-Top(B) TDC
  T->SetBranchAddress("R.s2.rt_c",Rs2rt_c); // Right arm S2-Bottom(A) TDC
 
  // Trigger condtion //
  T->SetBranchAddress("DR.evtypebits",&DR_evtypebits); 


////////////////////////////////////////////////////////////////////////


  //=======================================================//
  //=========== Set Parameters ============================//
  //======================================================//

 //========== Defined Parameters =======================//

  //---- TCanvas --------//
 TCanvas* cwo[chmax];
 TCanvas* cw[chmax];
 TCanvas* ctof[chmax];
 TCanvas* cital[chmax];
 //----- Function ------//
 TF1* ffit[chmax][2];
 TF1* ftof[chmax];
 TF1* ftof_c[chmax];
 TF1* ftof_min_c[chmax];
 TH2F* htof_adc[chmax][2];
 TH2F* htof_adc_c[chmax][2][ital1][ital2];
 TH2F* htdc_adc[chmax][2];
 TH2F* htdc_adc_c[chmax][2];
 TH1F* htdc[chmax][2];
 TH1F* htdc_c[chmax][2];
 TF1* ftdc[chmax][2];
 TF1* ftdc_c[chmax][2];
 TF1 *ftof_gaus[chmax];
 TH1F* htof[chmax];
 TH1F* htof_c[chmax][ital1][ital2];
 TH1F* htof_min_c[chmax];
 TH1F* htdc_m[chmax];
 TH1D* hmean_tof[chmax][2];
 TH2F* hital[chmax];
 //--- Range Parameters -----//
 double min_tof[chmax],max_tof[chmax],min_adc[chmax],max_adc[chmax];
 double bin_tof[chmax];
 double bin_adc[chmax];
 double bin_tdc[chmax];
 double p0[chmax][2],p1[chmax][2];
 double min_tof_c[chmax],max_tof_c[chmax],min_adc_c[chmax],max_adc_c[chmax];
 double bin_tof_c[chmax];
 double bin_adc_c[chmax];
 //---- Redefinition TBranch ------// 
 double adc[chmax][2];
 double tdc[chmax][2];
 double tof[chmax];
 double tdc_c[chmax][2];
 double tof_c[chmax];
 //----- Cut Definition ------//
 bool cut_time[chmax];
 bool cut_trig;
 bool cut_adc[chmax][2];
 //----- TOF Fit Parameters ---//
 double sig[chmax];
 double sig_c[chmax];
 double sig_min_c[chmax];
 for(int k=0;k<chmax;k++){
sig_c[k]=100.;
}
 double mean[chmax];
 double mean_c[chmax];
 //------ Itallation Weight ---//

double w0r,w0l,w2r[chmax],w2l[chmax];
 double sig_ital[chmax][ital1][ital2];
 double ws_min[chmax],wt_min[chmax];
 double wsmin[chmax],wtmin[chmax];

 int s_start=0; int s_end=ital1;
 int t_start=0; int t_end=ital2;
 int smin[chmax],tmin[chmax];


   //============ Parameter Setting =====================//

 for(int i=0;i<chmax;i++){  
       for(int j=0;j<2;j++){
 
 min_tof[i]=range_para(i,0);
 max_tof[i]=range_para(i,1);
 bin_tof[i]=(range_para(i,1)-range_para(i,0))/tdc_time;
 bin_tof[i]=(int)bin_tof[i];
 min_tof_c[i]=range_para(i,2);
 max_tof_c[i]=range_para(i,3);
 bin_tof_c[i]=(range_para(i,3)-range_para(i,2))/tdc_time;
 bin_tof_c[i]=(int)bin_tof_c[i];
 min_adc[i]=range_para(i,4);
 max_adc[i]=range_para(i,5); 
 bin_adc[i]=(-range_para(i,4)+range_para(i,5));
 bin_adc[i]=(int)bin_adc[i];
 min_adc_c[i]=range_para(i,4);
 max_adc_c[i]=range_para(i,5); 
 bin_adc_c[i]=(-range_para(i,4)+range_para(i,5));
 bin_adc_c[i]=(int)bin_adc_c[i];
       
  // Definition Function //
   if(j==0){     
  //------- TOF hist w/o Time-Walk Correction ---------//
  htof[i]=new TH1F(Form("htof[%d]",i),Form("TOF with Scinti CH %d, w/o correlation ",i),bin_tof[i],min_tof[i],max_tof[i]);
 //---------------- Fit TOF Hist --------------------//
ftof[i]=new TF1(Form("ftof[%d]",i),"gaus",min_tof[i],max_tof[i]);
ftof_c[i]=new TF1(Form("ftof_c[%d]",i),"gaus",min_tof_c[i],max_tof_c[i]);
ftof_min_c[i]=new TF1(Form("ftof_min_c[%d]",i),"gaus",min_tof_c[i],max_tof_c[i]);
 }
  //------------- TOF vs ADC FIt Function ------------//
 ffit[i][j]=new TF1(Form("ffit[%d][%d]",i,j),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc[i],max_adc[i]);
 //-------------- TOF Projection FitSlicesY ----------//
  ftof_gaus[i]=new TF1(Form("ftof_gaus[%d]",i),"gaus",min_tof[i],max_tof[i]);
  //---- TOF vs ADC 2D Hist w/o Time-Walk Correction ---//
 htof_adc[i][j]=new TH2F(Form("htof_adc[%d][%d]",i,j),Form("TOF vs Scinti CH %d,[%d]-PMT adc w/o correlation ",i,j),bin_adc[i],min_adc[i],max_adc[i],bin_tof[i],min_tof[i],max_tof[i]);

   //------------ TOF italtion  hist -------------//
   hital[i]=new TH2F(Form("hital[%d]",i),Form("TOF Itallation hist ch %d",i),ital1,0,ital1,ital2,0,ital2);

       }
        }           





  	
//========== Fill TOF vs ADC w/o Time-Walk CorrectionHist ==============//

 int ent;
 ent=T->GetEntries();
 cout<<"Event number : "<<ent<<endl;
 for(int k=0;k<ent;k++){ // Fill Entr
   T->GetEntry(k);

 for(int i=0;i<chmax;i++){  
   // if(chselect){
        if(i> ch_st-1 && i<ch_end+1){
   for(int j=0;j<2;j++){

    // ADC & TDC inf//
   adc[i][0]=Rs2ra_c[i];
   adc[i][1]=Rs2la_c[i];
   adc[16][0]=Rs0ra_c[0];
   adc[16][1]=Rs0la_c[0];
   tdc[i][0]=Rs2rt_c[i];
   tdc[i][1]=Rs2lt_c[i];
   tdc[16][0]=Rs0rt_c[0];
   tdc[16][1]=Rs0lt_c[0];


  //----- Cut Condition Difinition ------//
 cut_time[i]=false;
 cut_time[16]=false;
 cut_trig=false;
 if(DR_evtypebits==112)cut_trig=true;
 if(tdc[16][0]>0 && tdc[16][1]>0)cut_time[16]=true;
 if(tdc[i][0]>0 && tdc[i][1]>0)cut_time[i]=true;
 
 
  tof[i]=(tdc[i][0]+tdc[i][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0;

  //------- TOF w/0 Time-Walk Correction Hist --------//  
  if(j==0 && cut_trig && cut_time[i] && cut_time[16]){
  htof[i]->Fill(tof[i]);
  }
  //------- TOF vs ADC 2D Hist ----------------------//  
   if(cut_trig && cut_time[i] && cut_time[16])htof_adc[i][j]->Fill(adc[i][j],tof[i]); 
   
    
    }
      }
       }     
 }// END Fill


 //============================================================//
 //========== Time Walk correction ============================//
 //============================================================//


 //============ Get Time-Walk Correction Parameters ================//

 for(int i=0;i<chmax;i++){  
   // if(chselect){ 
     if(i> ch_st-1 && i<ch_end+1 || i==16){
   for(int j=0;j<2;j++){

 
 //------ FitSlicesY -----------------// 
 TObjArray aSlices;     
 htof_adc[i][j]->FitSlicesY(ftof_gaus[i],0,-1,0,"QRG3");
 hmean_tof[i][j]=(TH1D*)gROOT->FindObject(Form("htof_adc[%d][%d]_1",i,j));



 ffit[i][j]->SetParameter(0,fit_ini(i,j,0));   
 ffit[i][j]->SetParameter(1,fit_ini(i,j,1));   
 // ffit[i][j]->FixParameter(0,fit_ini(i,j,0));   
 // ffit[i][j]->FixParameter(1,1.0e8);   
 // ffit[i][j]->SetParLimits(0,0.5*fit_ini(i,j,0),1.5*fit_ini(i,j,0));   
 //ffit[i][j]->SetParLimits(1,0.5*fit_ini(i,j,1),1.5*fit_ini(i,j,1));   

 hmean_tof[i][j]->Fit(Form("ffit[%d][%d]",i,j),"","");

      p0[i][j]=ffit[i][j]->GetParameter(0);
      p1[i][j]=ffit[i][j]->GetParameter(1);
      //cout<<"p0["<<i<<"] : "<<p0[i][j]<<endl; 
      //cout<<"p1["<<i<<"] : "<<p1[i][j]<<endl; 

      //}
      }
 }     
 }

 
 //======== Fill TOF vs ADC 2D Hist w/ Time-Walk Correction ===========//


 //====== Func Definition ========//

 for(int s=s_start;s<s_end;s++){//Itallation 1
   for(int t=t_start;t<t_end;t++){ //Italltion 2
     for(int i=0;i<chmax;i++){  
      if(i> ch_st-1 && i<ch_end+1){
//---- TOF Hist w/ Time-Walk Correction ---//
 htof_c[i][s][t]=new TH1F(Form("htof_c[%d][%d][%d]",i,s,t),Form("TOF with Scinti CH %d, w/ correlation ",i),bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
 for(int j=0;j<2;j++){  
  //---- TOF vs ADC 2D Hist w/ Time-Walk Correction ---//
   htof_adc_c[i][j][s][t]=new TH2F(Form("htof_adc_c[%d][%d][%d][%d]",i,j,s,t),Form("TOF vs Scinti CH %d,[%d]-PMT adc w/ correlation ",i,j),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
 //}
	}
     }
   }
 }
 }
 // for(int s=s_start;s<s_end;s++){//Itallation 1
  //for(int t=t_start;t<t_end;t++){ //Italltion 2
    for(int i=0;i<chmax;i++){ 
        if(i> ch_st-1 && i<ch_end+1){

    for(int s=s_start;s<s_end;s++){//Itallation 1
  for(int t=t_start;t<t_end;t++){ //Italltion 2
        for(int k=0;k<ent;k++){ // Fill Entr
   T->GetEntry(k);
   // ADC & TDC inf//
   adc[i][0]=Rs2ra_c[i];
   adc[i][1]=Rs2la_c[i];
   adc[16][0]=Rs0ra_c[0];
   adc[16][1]=Rs0la_c[0];
   tdc[i][0]=Rs2rt_c[i];
   tdc[i][1]=Rs2lt_c[i];
   tdc[16][0]=Rs0rt_c[0];
   tdc[16][1]=Rs0lt_c[0];
    
     
  //----- Cut Condition Difinition ------//
 cut_time[i]=false;
 cut_time[16]=false;
 cut_trig=false;
 //cut_adc[i][j]=false;
 cut_adc[i][0]=false;
 cut_adc[i][1]=false;
 cut_adc[16][0]=false;
 cut_adc[16][1]=false;
 /*
 tdc_c[i][0]=0.0;  
 tdc_c[i][1]=0.0;
 tdc_c[16][0]=0.0;
 tdc_c[16][1]=0.0;
 */
 if(DR_evtypebits==112)cut_trig=true;
 if(tdc[16][0]>0 && tdc[16][1]>0)cut_time[16]=true;
 if(tdc[i][0]>0 && tdc[i][1]>0)cut_time[i]=true;
 
 if(adc[i][0]>100){
 cut_adc[i][0]=true;
 tdc_c[i][0]=(p0[i][0])*(1./sqrt(adc[i][0])-1./sqrt(p1[i][0]));
 }
 if(adc[i][1]>100){
 cut_adc[i][1]=true; 
tdc_c[i][1]=(p0[i][1])*(1./sqrt(adc[i][1])-1./sqrt(p1[i][1])); 
}
 if(adc[16][0]>1000){
 cut_adc[16][0]=true;
tdc_c[16][0]=(p0[16][0])*(1./sqrt(adc[16][0])-1./sqrt(p1[16][0]));
 }
 if(adc[16][1]>1000){
 cut_adc[16][1]=true;
 tdc_c[16][1]=(p0[16][1])*(1./sqrt(adc[16][1])-1./sqrt(p1[16][1]));
 }
 

   //-------- TOF w/ Time-Walk Correction Values -------------//
 //--------- TOF hist w/ Time-Walk Correction --------//
   // Itallation Weight Definition //
   w2r[i]=0.5+0.1*s;
   w2l[i]=0.5+0.1*t;
   w0r=0.888; //w0r 0.888 
   w0l=1.083; //w0l 1.083
   tof_c[i]=(tdc[i][0]+tdc[i][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0-w2r[i]*tdc_c[i][0]-w2l[i]*tdc_c[i][1]-w0r*tdc_c[16][0]-w0l*tdc_c[16][1];

  //------- TOF w/ Time-Walk Correction Hist --------//  
   if(cut_trig && cut_time[i] && cut_time[16] && cut_adc[i][0]&& cut_adc[i][1] && cut_adc[16][0] && cut_adc[16][1]){
  htof_c[i][s][t]->Fill(tof_c[i]);

   }
  //------- TOF w/ TW Corr vs ADC 2D Hist --------//  
       for(int j=0;j<2;j++){
   if(cut_trig && cut_time[i] && cut_time[16] && cut_adc[i][j]){
htof_adc_c[i][j][s][t]->Fill(adc[i][j],tof_c[i]); 
   }
       }
	}// END Fill Entr     
 
 //}
       
  

 htof_c[i][s][t]->Fit(ftof_c[i]);
 sig_ital[i][s][t]=ftof_c[i]->GetParameter(2);
 if(sig_c[i]>sig_ital[i][s][t]){
 sig_c[i]=sig_ital[i][s][t]; 
 smin[i]=s; tmin[i]=t;
 wsmin[i]=w2r[i];wtmin[i]=w2l[i];
 //htof_min_c[ch]=(TH1F*)htof_c[ch][s][t]->Clone();
 }
 cout<<"i "<<i<<": ital1 "<<s+1<<Form("/%d :ital2 ",ital1)<<t+1<<Form("/%d",ital2)<<endl;
 cout<<"sig_ital :"<<sig_ital[i][s][t]<<endl;
 cout<<"sig_c :"<<sig_c[i]<<endl;
 cout<<"smin "<<smin[i]<<" : tmin "<<tmin[i]<<endl;
 cout<<"wsmin "<<wsmin[i]<<" : wtmin "<<wtmin[i]<<endl;
 hital[i]->Fill(s,t,sig_ital[i][s][t]);
 int ital_entry=htof_c[i][s][t]->GetEntries();
 cout<<"Get Entries :"<<ital_entry<<endl;
       
  }
    }
	}
    }

//========= TOF Get Parameters ===========//
for(int i=0;i<chmax;i++){
        if(i> ch_st-1 && i<ch_end+1){     
 htof[i]->Fit(ftof[i]);
 mean[i]=ftof[i]->GetParameter(1);
 sig[i]=ftof[i]->GetParameter(2);

 htof_c[i][smin[i]][tmin[i]]->Fit(ftof_min_c[i]);
 mean_c[i]=ftof_min_c[i]->GetParameter(1);
 sig_min_c[i]=sig_c[i];
  }
 }
/////////////////////////////////////////


// int ital_ent=htof_c[i][smin][tmin]->GetEntries();
   //=ftof_min_c[ch]->GetParameter(2);

 //htof_c[ch]->Fit(ftof_c[ch]);
 // mean_c[ch]=ftof_c[ch]->GetParameter(1);
 // sig_c[ch]=ftof_c[ch]->GetParameter(2);


 //================================================//
 //======== Draw Hist ============================//
 //===============================================//

 bool print=true;
 bool write=true;
 bool draw=true;
 bool comment=true;
 
 int print_i=0;
 
 for(int i=0;i<chmax;i++){   
  int ch=i;
          if(i> ch_st-1 && i<ch_end+1){
// Definition TCanvas //
 
     if(draw ||print){
       cwo[ch]=new TCanvas(Form("cwo[%d]",ch),Form("cwo[%d]",ch));
       cwo[ch]->Divide(1,2);
       cw[ch]=new TCanvas(Form("cw[%d]",ch),Form("cw[%d]",ch));
       cw[ch]->Divide(1,2);
       ctof[ch]=new TCanvas(Form("ctof[%d]",ch),Form("ctof[%d]",ch));
       ctof[ch]->Divide(1,2);
     }
       if(draw){
   for(int j=0;j<2;j++){
     // TOF vs ADC w/o Time-Walk Correction //
     cwo[ch]->cd(j+1);
     htof_adc[ch][j]->Draw();
     hmean_tof[ch][j]->Draw("same");
     // TOF vs ADC w/ Time-Walk Correction //
     cw[ch]->cd(j+1);
     htof_adc_c[ch][j][smin[ch]][tmin[ch]]->Draw();
       }     
 
 //  TOF hist w/o Time-Walk Correction//
     ctof[ch]->cd(1);
     htof[ch]->GetXaxis()->SetRangeUser(mean[ch]-5*sig[ch],mean[ch]+5*sig[ch]);
     htof[ch]->Draw();
     //  TOF hist w/ Time-Walk Correction//
     ctof[ch]->cd(2);
     htof_c[ch][smin[ch]][tmin[ch]]->Draw();
 //========= TOF Itallation Hist =========//
   cital[ch]=new TCanvas(Form("cital[%d]",ch),Form("cital[%d]",ch));
   cital[ch]->cd(1);
   hital[ch]->Draw("colz");
   }


//========= Comment Out ========================//
   if(comment){
 cout<<"======== Comment Out =============="<<endl;
 cout<<"TOF S2 CH"<<ch<<" w/o Time-Walk sigma : "<<sig[ch]<<endl;
 cout<<"TOF S2 CH"<<ch<<" w/ Time-Walk sigma : "<<sig_c[ch]<<endl;
 cout<<"p0[ch][0] : "<<p0[ch][0]<<" p0[ch][1] : "<<p0[ch][1]<<endl;
 cout<<"p0[16][0] : "<<p0[16][0]<<" p0[16][1] : "<<p0[16][1]<<endl;
 cout<<"p1[ch][0] : "<<p1[ch][0]<<" p1[ch][1] : "<<p1[ch][1]<<endl;
 cout<<"p1[16][0] : "<<p1[16][0]<<" p1[16][1] : "<<p1[16][1]<<endl;
 cout<<"smin "<<smin[ch]<<" : tmin "<<tmin[ch]<<endl;
 cout<<Form("sig_min_c [w2r=%lf][w02l=%lf] : ",wsmin[ch],wtmin[ch])<<sig_min_c[ch]<<endl;
 // cout<<"htof_min_c[ch] Event number : "<<ital_ent<<endl;

   }

 //===================== Write Out =================================//
 
 if(write){
 const char* fileName="./param/twcorr.param"; 
 std::ofstream ofs(fileName,std::ios::app); 

 char *buff;
 time_t now =std::time(NULL);
 struct tm* pnow=localtime(&now);
 buff=Form("%04d/%02d/%02d %02d:%02d, %02d sec",pnow->tm_year+1900,pnow->tm_mon+1,pnow->tm_mday,pnow->tm_hour,pnow->tm_min,pnow->tm_sec);
 
if(!ofs){
   std::cout<<"Couldn't open file "<<std::endl;
   std::cin.get();
   return 0;
} 

 ofs<<"===== <"<<buff<<" > ==========="<<std::endl;
 ofs<<" CH "<<" PMT " <<" p0 "<<" p1 "<<" S2 weight "<<" TOF sig "<<std::endl;
 ofs<<" "<< ch <<" "<< 0 <<" "<< p0[ch][0]<<" "<< p1[ch][0]<<" "<< wsmin[ch] <<" "<< sig_min_c[ch]<<std::endl;
ofs<<" "<< ch <<" "<< 1 <<" "<< p0[ch][1]<<" "<< p1[ch][1]<<" "<< wtmin[ch] <<" "<< sig_min_c[ch]<<std::endl;
 ofs<<"============== END ==================="<<std::endl;
 std::cin.get();
 }
 
	  }
 }



 //================ Print Canvas =================================//

TString name;
 if(print){
    for(int ch=ch_st;ch<ch_end+1;ch++){
   if(ch==ch_st){
 name.Form("./pdf/twcorr_pdf/tof_twcorr.pdf");
 ctof[ch]->Print(name+"[","pdf");
   }
 ctof[ch]->Print(name,"pdf");
 cw[ch]->Print(name,"pdf");
 cwo[ch]->Print(name,"pdf");
 cital[ch]->Print(name,"pdf");
 if(ch==ch_end)cital[ch]->Print(name +"]","pdf");
}
 }
 

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


