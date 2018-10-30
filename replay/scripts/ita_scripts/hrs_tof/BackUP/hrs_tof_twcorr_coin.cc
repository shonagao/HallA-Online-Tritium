// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2


// Grobal Function //
 int chmax=17; // channel of S2 PMT 
 int ch=8;// defolt ch 8;
 double tdc_time=0.1e-9;//TDC converse ch->sec [ss/ch]
// Itallation parameters
int ital1=10;
int ital2=10;
  Double_t tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec

//====== Define Function ===========//
double range_para(int i,int j);
double range_adc(int i,int j);
double fit_ini(int i,int j,int k);

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
/*
int main(int argc, char** argv){
  TApplication *theApp =new TApplication("App",&argc,argv);
*/

void hrs_tof_twcorr_coin(){

//-------- TTree data input ---------------//
  TChain*  T=new TChain("T");
    for(Int_t i = 93498; i<93500;i++) /// coincidence run goes 100400 to 100610
  //  for(Int_t i = 100400; i<100610;i++) /// coincidence run goes 100400 to 100610
 {
   //   T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
   T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d*.root",i));//For coincidence run
 }

 //=================================================//
 //========== Get Tree Branch ======================//
 //=================================================//

 // Variable to be used event by event 
	  Double_t RF_s2_mean;
   // Define Tree leaf variable to hold the values
            TString nHRS_trig = "DR.rrawt4";
	  //TString nHRS_trig = "DR.evtypebits"; // for beta =1 particle
          TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
	  TString nRs2_tdchit = "RTDC.F1FirstHit"; // for coincidence run
	  //  TString nRs2_tdchit = "F1FirstHit"; // for marathan run
	  TString nR_trx = "R.tr.x";
          TString nR_trth = "R.tr.th";
          TString nRs2_lac = "R.s2.la_c";
          TString nR_trbeta="R.tr.beta";

    // Define Tree leaf variable
	  Double_t HRS_trig;
	  Double_t Rs2_pads[100];
	  Double_t Rs2_nthit;
	  Double_t Rs2_tdchit[100];
	  Double_t R_trx[100];
	  Double_t R_trth[100];
	  Double_t Rs2_lac[100];
          Double_t R_trbeta;
    // Set The Branch addres

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
  T->SetBranchStatus(nHRS_trig,1);
  T->SetBranchStatus(nRs2_pads,1);
  T->SetBranchStatus(nRs2_nthit,1);
  T->SetBranchStatus(nRs2_tdchit,1);
  T->SetBranchStatus(nR_trx,1);
  T->SetBranchStatus(nR_trth,1);
  T->SetBranchStatus(nRs2_lac,1);
  T->SetBranchStatus(nR_trbeta,1);

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

  //=== Coincidence F1TDC Branch ================//
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);
  //============================================//

////////////////////////////////////////////////////////////////////////


	  // original values
 Double_t corr_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};
 /// the following loop changes from Marathon to coincidence runs
 /*
 Double_t c =  -3.9094e-09;  // -3.01594e-9 // real value -3.9094e-09; 
 for(int i=0;i<14;i++)
   {
     corr_x[i] = corr_x[i] + c;
   }
 */
 // since I am using the coincidence data now and there is some offset. I want to remove the offset. So I would like to add a constant in each of the array element.	 
 Double_t corr_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};	 
 Double_t corr_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};
 // the following corrections will not be applied

 // Double_t alignment[14] = {2.0975e-9, 2.1e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
   Double_t alignment[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
   // This value is taken from coin_trig.C which is the common trig for LHRS and RHrs. Actually I took this value from test_coin.C


        /// loop over entries
/////////////////////////////////////////////////////////////

	      /// making loop for each paddle
	      /* 




	  Long64_t nentries = T->GetEntries();
	  for(Long64_t j=0;j<nentries;j++)
	    {
	      T->GetEntry(j); 
	    }


for( int i=2;i<16;i++)
    {
      if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01  ) // R.tr.beta>0.76 && R.tr.beta<0.80
	{
	  //Calc RF_S2mean depending on i
	  double RF_s2_mean = (Rs2_tdchit[8] - Rs2_tdchit[1]) * tdcTime  - (((Rs2_tdchit[8]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0 + corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2]);
	    


	   h1->Fill(RF_s2_mean);// for 2 Dimensional;
	  //  h1->Fill(RF_s2_mean,i);// For opne  dimensional

 
	}
      // else if ( HRS_trig>1 && Rs2_nthit==2 && Rs2_pads[0]==i && Rs2_pads[1]==i+1 &&  R_trbeta>0.76 && R_trbeta<80) //  R.tr.beta>0.76 && R.tr.beta<80
      //	{
      //	}
      // else {}
    }


////////////////////////////////////////////////////////

} // end entry loop ie irun=90854 loop

    }// end subrun
    }

  // TF1 *f1 = new TF1("f1","[0]+[1]*x",-5,12);
  // h1->Fit("f1","R+");
 
  TCanvas *c1 = new TCanvas("c1", "c1",600,600);
  c1->cd();
  h1->Draw("colz");
  // p1->Add(h1);
  // For title and some other headings
 
  h1->SetTitle("R-HRS S2 Timing Alignments");
  gStyle->SetTitleY(0.96);
 
  h1->GetXaxis()->SetTitle("RF-S2 Meantime(sec)");
  h1->GetXaxis()->SetTitleSize(0.035);
  h1->GetXaxis()->SetTitleOffset(1.02);
  h1->GetXaxis()->CenterTitle();
  
  h1->GetYaxis()->SetTitle("Counts");
  h1->GetYaxis()->CenterTitle();
  h1->GetYaxis()->SetTitleOffset(1.02);
  
  TLatex l1;
  l1.SetTextSize(0.032);
  l1.DrawLatex(0.575e-6,700,Form("#color[2]{#sigma #approx 395 ps}"));
 
   TLine *l2= new TLine(0.325e-6,0.5,0.325e-6,14.5);
   l2->SetLineColor(kRed);
   l2->Draw();
  //c1->SaveAs("RHRStime.pdf"); 


  */





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
 bool cut_hits[chmax];
bool cut_coin[chmax];
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
 int smin,tmin;

//--------- Coincidence run ------------//
double rftdc;
double corr[chmax];


   //============ Parameter Setting =====================//
 for(int i=0;i<chmax;i++){  
      if(i==ch){
   for(int j=0;j<2;j++){
 
 min_tof[i]=range_para(i,0);
 max_tof[i]=range_para(i,1);

 min_tof[i]=-1.0e-5;
 max_tof[i]=1.0e-5;
 bin_tof[i]=(range_para(i,1)-range_para(i,0))/tdcTime;
 
 bin_tof[i]=(int)bin_tof[i];
 min_tof_c[i]=range_para(i,2);
 max_tof_c[i]=range_para(i,3);
 bin_tof_c[i]=(range_para(i,3)-range_para(i,2))/tdcTime;
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
      if(i==ch){
   for(int j=0;j<2;j++){

    // ADC & TDC inf Single ARM condition//
     /*
   adc[i][0]=Rs2ra_c[i];
   adc[i][1]=Rs2la_c[i];
   adc[16][0]=Rs0ra_c[0];
   adc[16][1]=Rs0la_c[0];
   tdc[i][0]=Rs2rt_c[i];
   tdc[i][1]=Rs2lt_c[i];
   tdc[16][0]=Rs0rt_c[0];
   tdc[16][1]=Rs0lt_c[0];
     */

    // ADC & TDC inf Coincidence condition//
   adc[i][0]=Rs2ra_c[i];
   adc[i][1]=Rs2la_c[i];
   tdc[i][0]=(Rs2_tdchit[8]-Rs2_tdchit[i+16])*tdcTime;
   tdc[i][1]=(Rs2_tdchit[48]-Rs2_tdchit[i+48])*tdcTime;
   rftdc=(Rs2_tdchit[8]-Rs2_tdchit[1])*tdcTime;
corr[i]=corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2];


 tof[i]=rftdc-(tdc[i][0]+tdc[i][1])/2.0+corr[i]; //Coin

  //----- Cut Condition Difinition ------//
 cut_time[ch]=false;
 cut_time[16]=false;
 cut_trig=false;
 cut_coin[i]=false;
 if(DR_evtypebits==112)cut_trig=true;
 if(tdc[16][0]>0 && tdc[16][1]>0)cut_time[16]=true;
 if(tdc[ch][0]>0 && tdc[ch][1]>0)cut_time[ch]=true;
 if(HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01 )cut_coin[i]=true;

  //------- TOF w/0 Time-Walk Correction Hist --------//  
 if(HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01 && j==0 ){
  htof[i]->Fill(tof[i]);
  }
  //------- TOF vs ADC 2D Hist ----------------------//  
    if(HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01 )htof_adc[i][j]->Fill(adc[i][j],tof[i]); 
   

     }
      }
 }     

 
 }// END Fill

 
 for(int i=0;i<chmax;i++){  
      if(i==ch){
   for(int j=0;j<2;j++){

     // TOF vs ADC w/o Time-Walk Correction //
     cwo[i]->cd(j+1);
     htof_adc[i][j]->Draw();
     //     hmean_tof[i][j]->Draw("same");
  
     }
      }
       }     
 


     ctof[ch]=new TCanvas(Form("ctof[%d]",ch),Form("ctof[%d]",ch));
     //ctof[ch]->Divide(1,2);
     ctof[ch]->cd(1);
     // htof[ch]->GetXaxis()->SetRangeUser(mean[ch]-5*sig[ch],mean[ch]+5*sig[ch]);
     htof[ch]->Draw();
/*

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
 
 if(adc[ch][0]>100){
 cut_adc[ch][0]=true;
 tdc_c[ch][0]=(p0[ch][0])*(1./sqrt(adc[ch][0])-1./sqrt(p1[ch][0]));
 }
 if(adc[ch][1]>100){
 cut_adc[ch][1]=true; 
tdc_c[ch][1]=(p0[ch][1])*(1./sqrt(adc[ch][1])-1./sqrt(p1[ch][1])); 
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
   w2r[ch]=0.0+0.1*s;
   w2l[ch]=0.0+0.1*t;
   w0r=0.888; //w0r 0.888 
   w0l=1.083; //w0l 1.083
   tof_c[i]=(tdc[ch][0]+tdc[ch][1])/2.0-(tdc[16][0]+tdc[16][1])/2.0-w2r[ch]*tdc_c[ch][0]-w2l[ch]*tdc_c[ch][1]-w0r*tdc_c[16][0]-w0l*tdc_c[16][1];

  //------- TOF w/ Time-Walk Correction Hist --------//  
   if(i==ch && j==0 &&cut_trig && cut_time[ch] && cut_time[16] && cut_adc[ch][0]&& cut_adc[ch][1] && cut_adc[16][0] && cut_adc[16][1]){
  htof_c[ch][s][t]->Fill(tof_c[i]);
   }
  //------- TOF w/ TW Corr vs ADC 2D Hist --------//  
   if(cut_trig && cut_time[ch] && cut_time[16] && cut_adc[i][j]){
htof_adc_c[i][j][s][t]->Fill(adc[i][j],tof_c[i]); 
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
 wsmin=w2r[ch];wtmin=w2l[ch];
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
 cout<<"smin "<<smin<<" : tmin "<<tmin<<endl;
 cout<<Form("sig_min_c [w2r=%lf][w02l=%lf] : ",wsmin,wtmin)<<sig_min_c[ch]<<endl;
 cout<<"htof_min_c[ch] Event number : "<<ital_ent<<endl;

*/

//theApp->Run();
// return 0;



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






 double fit_ini(int i,int j,int k){
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


