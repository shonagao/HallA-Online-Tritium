// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2


// Grobal Function //
 int chmax=17; // channel of S2 PMT 
 int ch=8;// defolt ch 8;
 double tdc_time=0.1e-9;//TDC converse ch->sec [ss/ch]
// Itallation parameters
int ital1=10;
int ital2=10;


//====== Define Function ===========//
double range_para(int i,int j);
double range_adc(int i,int j);
long double fit_ini(int i,int j,int k);

//=================================//
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
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
 for(Int_t i = 93500; i<93501;i++) /// these run goes 93495 to 93504
 {
    T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }

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
  T->SetBranchAddress("R.s2.lt_c",Rs2lt_c); // Right arm S2-Top(B) TDC
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
 sig_c[ch]=100.;
 double mean[chmax];
 double mean_c[chmax];
 //------ Itallation Weight ---//

double w0r,w0l,w2r[chmax],w2l[chmax];
 double sig_ital[chmax][ital1][ital2];
 double ws_min,wt_min;
 double wsmin,wtmin;

 int s_start=0; int s_end=ital1;
 int t_start=0; int t_end=ital2;

// Itallation Parameters Setting //

/*
//Weight Definition ( weight= wscen+ ws*s)
 double ws_start=0.0;
 double wt_start=0.0; //center of weight value examp.  
 double ws_end=2.0;
 double wt_end=2.0; 

 
 double ws=(ws_end-ws_start)/ital1_d;
 double ws=(wt_end-wt_start)/ital2_d;
double ital1=ital1_d;
double ital2=ital2_d;
ital1=(int)ital1;// itallation counts
ital2=(int)ital2;

 // double wt=0.1; // weight intarval
 int s_start,s_end,t_start,t_end;
 s_start=0; s_end=ital1;
 t_start=0; t_end=ital2;
 double wsmin,wsmax,wtmin,wtmax;
 wsmin=ws_start;
 wsmax=ws_end;
 wtmin=wt_start;
 wtmax=wt_end;
 double wsbin,wtbin;
 wsbin=(wsmax-wsmin);
 wsbin=(int)wsbin;
 wtbin=(wtmax-wtmin);
 wtbin=(int)wtbin;

*/ 
 //------ Chi squire ---------//

 // double chi2[chmax][ital1][ital2];
 //double chi2_min[chmax];
 // chi2_min[ch]=10000.0;
 int smin,tmin;
   //============ Parameter Setting =====================//
 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){
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


 if(i==ch && j==0){     
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

      }
   //------------ TOF italtion  hist -------------//
   hital[i]=new TH2F(Form("hital[%d]",i),Form("TOF Itallation hist ch %d",i),ital1,0,ital1,ital2,0,ital2);
// Definition TCanvas //
       cwo[i]=new TCanvas(Form("cwo[%d]",i),Form("cwo[%d]",i));
                         cwo[i]->Divide(1,2);
       cw[i]=new TCanvas(Form("cw[%d]",i),Form("cw[%d]",i));
                        
      }
 }

 	
//========== Fill TOF vs ADC w/o Time-Walk CorrectionHist ==============//

 int ent;
 ent=T->GetEntries();
 cout<<"Event number : "<<ent<<endl;
 for(int k=0;k<ent;k++){ // Fill Entr
   T->GetEntry(k);

 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){
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
 cut_time[ch]=false;
 cut_time[16]=false;
 cut_trig=false;
 if(DR_evtypebits==112)cut_trig=true;
 if(tdc[16][0]>0 && tdc[16][1]>0)cut_time[16]=true;
 if(tdc[ch][0]>0 && tdc[ch][1]>0)cut_time[ch]=true;
 
 
  tof[i]=(tdc[ch][0]+tdc[ch][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0;

  //------- TOF w/0 Time-Walk Correction Hist --------//  
  if(i==ch && j==0 && cut_trig && cut_time[ch] && cut_time[16]){
  htof[ch]->Fill(tof[ch]);
  }
  //------- TOF vs ADC 2D Hist ----------------------//  
   if(cut_trig && cut_time[ch] && cut_time[16])htof_adc[i][j]->Fill(adc[i][j],tof[i]); 
   






     }
      }
       }     

 


         }// END Fill




 //============================================================//
 //========== Time Walk correction ============================//
 //============================================================//


 //============ Get Time-Walk Correction Parameters ================//

 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){
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

    
     }
      }
       }     


 
 //======== Fill TOF vs ADC 2D Hist w/ Time-Walk Correction ===========//


 //====== Func Definition ========//

 for(int s=s_start;s<s_end;s++){//Itallation 1
   for(int t=t_start;t<t_end;t++){ //Italltion 2
     for(int i=0;i<chmax;i++){  
        if(i==ch || i==16){
//---- TOF Hist w/ Time-Walk Correction ---//
 htof_c[i][s][t]=new TH1F(Form("htof_c[%d][%d][%d]",i,s,t),Form("TOF with Scinti CH %d, w/ correlation ",i),bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
 for(int j=0;j<2;j++){  
  //---- TOF vs ADC 2D Hist w/ Time-Walk Correction ---//
 htof_adc_c[i][j][s][t]=new TH2F(Form("htof_adc_c[%d][%d][%d][%d]",i,j,s,t),Form("TOF vs Scinti CH %d,[%d]-PMT adc w/ correlation ",i,j),bin_adc_c[i],min_adc_c[i],max_adc_c[i],bin_tof_c[i],min_tof_c[i],max_tof_c[i]);
         }
	}
     }
   }
 }

 for(int s=s_start;s<s_end;s++){//Itallation 1
   for(int t=t_start;t<t_end;t++){ //Italltion 2
 for(int k=0;k<ent;k++){ // Fill Entr
   T->GetEntry(k);
 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){  
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
 cut_time[ch]=false;
 cut_time[16]=false;
 cut_trig=false;
 //cut_adc[i][j]=false;
 cut_adc[ch][0]=false;
 cut_adc[ch][1]=false;
 cut_adc[16][0]=false;
 cut_adc[16][1]=false;
 tdc_c[ch][0]=0.0;  
 tdc_c[ch][1]=0.0;
 tdc_c[16][0]=0.0;
 tdc_c[16][1]=0.0;
 if(DR_evtypebits==112)cut_trig=true;
 if(tdc[16][0]>0 && tdc[16][1]>0)cut_time[16]=true;
 if(tdc[ch][0]>0 && tdc[ch][1]>0)cut_time[ch]=true;
 
 if(adc[ch][0]>300){
 cut_adc[ch][0]=true;
 tdc_c[ch][0]=(p0[ch][0])*(1./sqrt(adc[ch][0])-1./sqrt(p1[ch][0]));
 }
 if(adc[ch][1]>300){
 cut_adc[ch][1]=true; 
tdc_c[ch][1]=(p0[ch][1])*(1./sqrt(adc[ch][1])-1./sqrt(p1[ch][1])); 
}
 if(adc[16][0]>100){
 cut_adc[16][0]=true;
tdc_c[16][0]=(p0[16][0])*(1./sqrt(adc[16][0])-1./sqrt(p1[16][0]));
 }
 if(adc[16][1]>100){
 cut_adc[16][1]=true;
 tdc_c[16][1]=(p0[16][1])*(1./sqrt(adc[16][1])-1./sqrt(p1[16][1]));
 }
 

   //-------- TOF w/ Time-Walk Correction Values -------------//
 //--------- TOF hist w/ Time-Walk Correction --------//
   // Itallation Weight Definition //
   w2r[ch]=0.5;
   w2l[ch]=0.5;
   w0r=0.88+0.001*s; //w0r 0.886 
   w0l=1.08+0.001*t; //w0l 1.087
   tof_c[i]=(tdc[ch][0]+tdc[ch][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0-w2r[ch]*tdc_c[ch][0]-w2l[ch]*tdc_c[ch][1]-w0r*tdc_c[16][0]-w0l*tdc_c[16][1];

   //     tof_c[i]=(tdc[ch][0]+tdc[ch][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0-w2r[ch]*tdc_c[ch][0]-w2l[ch]*tdc_c[ch][1]-w0r*tdc_c[16][0]-w0l*tdc_c[16][1];

  //------- TOF w/ Time-Walk Correction Hist --------//  
   if(i==ch && j==0 &&cut_trig && cut_time[ch] && cut_time[16] && cut_adc[ch][0]&& cut_adc[ch][1] && cut_adc[16][0] && cut_adc[16][1]){
  htof_c[ch][s][t]->Fill(tof_c[i]);
  
 if(k<100){
   cout<<"=========== Test ================"<<endl;
   cout<<"Cut conditon :"<<cut_trig<<":"<<cut_time[ch]<<":"<<cut_time[16]<<":"<<cut_adc[i][j]<<endl;
   cout<<"k is "<<k<<endl;
   cout<<"htof_c hist  (tof_c): "<<tof_c[i]<<endl;
   cout<<" i "<<i<<":j "<<j<<":s "<<s<<":t "<<t<<endl; 
}
   }

  //------- TOF w/ TW Corr vs ADC 2D Hist --------//  
   if(cut_trig && cut_time[ch] && cut_time[16] && cut_adc[i][j]){
     // tof_c[i]=(tdc[ch][0]+tdc[ch][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0-w2r[ch]*tdc_c[ch][0]-w2l[ch]*tdc_c[ch][1]-w0r*tdc_c[16][0]-w0l*tdc_c[16][1];
htof_adc_c[i][j][s][t]->Fill(adc[i][j],tof_c[i]); 
  
 if(k<100){
   cout<<"=========== Test ================"<<endl;
   cout<<"k is "<<k<<endl;
   cout<<"Cut conditon :"<<cut_trig<<":"<<cut_time[ch]<<":"<<cut_time[16]<<":"<<cut_adc[i][j]<<endl;
   cout<<"htof_adc_c : "<<tof_c[i]<<endl;
   cout<<" i="<<i<<" :j="<<j<<" :s="<<s<<" :t="<<t<<endl; 
}

   }
 }// END Fill Entr     

   }
      }
 }



 htof_c[ch][s][t]->Fit(ftof_c[ch]);
sig_ital[ch][s][t]=ftof_c[ch]->GetParameter(2);

 if(sig_c[ch]>sig_ital[ch][s][t]){
 sig_c[ch]=sig_ital[ch][s][t]; 
 smin=s; tmin=t;
 wsmin=w0r;wtmin=w0l;
 //htof_min_c[ch]=(TH1F*)htof_c[ch][s][t]->Clone();
 }
 cout<<"ital1 "<<s+1<<Form("/%d :ital2 ",ital1)<<t+1<<Form("/%d",ital2)<<endl;
 cout<<"sig_ital :"<<sig_ital[ch][s][t]<<endl;
 cout<<"sig_c :"<<sig_c[ch]<<endl;
 cout<<"smin "<<smin<<" : tmin "<<tmin<<endl;
 cout<<"wsmin "<<wsmin<<" : wtmin "<<wtmin<<endl;
 hital[ch]->Fill(s,t,sig_ital[ch][s][t]);
 int ital_entry=htof_c[ch][s][t]->GetEntries();
 cout<<"Get Entries :"<<ital_entry<<endl;

   }
 }



//========= TOF Get Parameters ===========//
 htof[ch]->Fit(ftof[ch]);
 mean[ch]=ftof[ch]->GetParameter(1);
 sig[ch]=ftof[ch]->GetParameter(2);
 htof_c[ch][smin][tmin]->Fit(ftof_min_c[ch]);
 mean_c[ch]=ftof_min_c[ch]->GetParameter(1);
 sig_min_c[ch]=sig_c[ch];
 int ital_ent=htof_c[ch][smin][tmin]->GetEntries();
   //=ftof_min_c[ch]->GetParameter(2);

 //htof_c[ch]->Fit(ftof_c[ch]);
 // mean_c[ch]=ftof_c[ch]->GetParameter(1);
 // sig_c[ch]=ftof_c[ch]->GetParameter(2);


 //================================================//
 //======== Draw Hist ============================//
 //===============================================//

 
 
 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){
   for(int j=0;j<2;j++){

     // TOF vs ADC w/o Time-Walk Correction //
     cwo[i]->cd(j+1);
     htof_adc[i][j]->Draw();
     hmean_tof[i][j]->Draw("same");
  
     }
      }
       }     
 

 
 for(int i=0;i<chmax;i++){  
      if(i==ch || i==16){
   cw[i]->Divide(1,2);
   for(int j=0;j<2;j++){
   // TOF vs ADC w/ Time-Walk Correction //

      cw[i]->cd(j+1);
      htof_adc_c[i][j][smin][tmin]->Draw();
     }
      }
       }     


 //  TOF hist w/o Time-Walk Correction//

     ctof[ch]=new TCanvas(Form("ctof[%d]",ch),Form("ctof[%d]",ch));
     ctof[ch]->Divide(1,2);
     ctof[ch]->cd(1);
     htof[ch]->GetXaxis()->SetRangeUser(mean[ch]-5*sig[ch],mean[ch]+5*sig[ch]);
     htof[ch]->Draw();
     //  TOF hist w/ Time-Walk Correction//
     ctof[ch]->cd(2);
     // htof_c[ch]->GetXaxis()->SetRangeUser(mean_c[ch]-5*sig_c[ch],mean_c[ch]+5*sig_c[ch]);
     //  htof_c[ch]->Draw();
     // htof_min_c[ch]->GetXaxis()->SetRangeUser(mean_c[ch]-5*sig_min_c[ch],mean_c[ch]+5*sig_min_c[ch]);
  htof_c[ch][smin][tmin]->Draw();
 //========= TOF Itallation Hist =========//
   
   cital[ch]=new TCanvas(Form("cital[%d]",ch),Form("cital[%d]",ch));
   cital[ch]->cd(1);
   hital[ch]->Draw("colz");
 



//==============================================//
//========= Comment Out ========================//
//==============================================//

 cout<<"======== Comment Out =============="<<endl;
 cout<<"TOF S2 CH"<<ch<<" w/o Time-Walk sigma : "<<sig[ch]<<endl;
 cout<<"TOF S2 CH"<<ch<<" w/ Time-Walk sigma : "<<sig_c[ch]<<endl;
 cout<<"p0[ch][0] : "<<p0[ch][0]<<" p0[ch][1] : "<<p0[ch][1]<<endl;
 cout<<"p0[16][0] : "<<p0[16][0]<<" p0[16][1] : "<<p0[16][1]<<endl;
 cout<<"p1[ch][0] : "<<p1[ch][0]<<" p1[ch][1] : "<<p1[ch][1]<<endl;
 cout<<"p1[16][0] : "<<p1[16][0]<<" p1[16][1] : "<<p1[16][1]<<endl;
 cout<<Form("sig_min_c [w0r=%lf][w0l=%lf] : ",wsmin,wtmin)<<sig_min_c[ch]<<endl;
 cout<<"htof_min_c[ch] Event number : "<<ital_ent<<endl;


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


