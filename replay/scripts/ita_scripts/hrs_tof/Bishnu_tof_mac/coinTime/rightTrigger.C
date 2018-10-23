//////////////////////////////////////////////////////////
//07/01/2018
//Author Bishnu Pandey
//Here I am going to combine all of the scintillator together except 0 & 1
//First two scintillator doesnot have enough statistics, so they are not included in this code
///////////////////////////////////////////////////////////////////


#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
 
Int_t NEvent = 1000;

void rightTrigger()
{
  Double_t tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec
  TH1F *h1 = new TH1F("h1", "RHRS s2 Detector common time(i<16);s2 time in sec;counts",500,0.55e-6,0.6e-6);
  TString filename;
  for(int irun= 90854;irun<90856;irun++)
    {

      if (irun == 90857) continue;


      for(int subrun =1; subrun<9; subrun++)
	{
	 
	  if( irun==90859 || subrun>>5) continue;

 filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d_%d.root", irun, subrun));

	   TFile *data_file = new TFile(filename,"READ");
           TTree *T = (TTree*)data_file->Get("T");
	 
          // Variable to be used event by event 
	  Double_t RF_s2_mean;
	  Double_t RF_s2_meanCorr_pad_2;
          Double_t RF_s2_meanCorr_pad_3;
          Double_t RF_s2_meanCorr_pad_4;
          Double_t RF_s2_meanCorr_pad_5;
          Double_t RF_s2_meanCorr_pad_6;
          Double_t RF_s2_meanCorr_pad_7;
          Double_t RF_s2_meanCorr_pad_8;
          Double_t RF_s2_meanCorr_pad_9;
          Double_t RF_s2_meanCorr_pad_10;
          Double_t RF_s2_meanCorr_pad_11;
          Double_t RF_s2_meanCorr_pad_12;
          Double_t RF_s2_meanCorr_pad_13;
          Double_t RF_s2_meanCorr_pad_14;
          Double_t RF_s2_meanCorr_pad_15;

	  // Define Tree leaf variable to hold the values
	  TString nHRS_trig = "DR.evtypebits";
          TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
          TString nRs2_tdchit = "F1FirstHit";
	  TString nR_trx = "R.tr.x";
          TString nR_trth = "R.tr.th";
          TString nRs2_lac = "R.s2.la_c";

	  // Define Tree leaf variable
	  Double_t HRS_trig;
	  Double_t Rs2_pads[100];
	  Double_t Rs2_nthit;
	  Double_t Rs2_tdchit[100];
	  Double_t R_trx[100];
	  Double_t R_trth[100];
	  Double_t Rs2_lac[100];

	  // Set The Branch addres
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);

	 
	  Double_t corr_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};

	  Double_t corr_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};

	  Double_t corr_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};

	  Double_t alignment[14] = {2.0975e-9, 2.1e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
 /// loop over entries

	  Long64_t nentries = T->GetEntries();
	  for(Long64_t i=0;i<nentries;i++)
	    {
	      T->GetEntry(i); 
	      UInt_t HRS_trig_bit5 =(static_cast<UInt_t>(HRS_trig)>>5)&1;
	      if( HRS_trig_bit5==0) continue;


 RF_s2_meanCorr_pad_2 = (Rs2_tdchit[15]-Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[18])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[50])*tdcTime)/2.0 + 4.91939e-09*R_trx[0] -4.46911e-08*R_trth[0] -8.66369e-13*Rs2_lac[2] + 2.0975e-9 );

	     
RF_s2_meanCorr_pad_3 = (Rs2_tdchit[15] - Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[19])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[51])*tdcTime)/2.0 +5.41034e-09 *R_trx[0] -4.79507e-08*R_trth[0] - 3.84824e-13*Rs2_lac[3] + 2.1e-9); 


RF_s2_meanCorr_pad_4 = (Rs2_tdchit[15]-Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[20])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[52])*tdcTime)/2.0 + 6.87688e-09*R_trx[0] -3.59540e-08*R_trth[0] -1.45016e-12*Rs2_lac[4] +0.85e-9 );


RF_s2_meanCorr_pad_5 = (Rs2_tdchit[15] - Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[21])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[53])*tdcTime)/2.0 + 9.22121e-09 *R_trx[0] -3.04303e-08*R_trth[0] - 12.08217e-13*Rs2_lac[5] +2.0e-9); 
      
RF_s2_meanCorr_pad_6 = (Rs2_tdchit[15]-Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[22])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[54])*tdcTime)/2.0  + 0.795016e-8*R_trx[0] - 3.25733e-08*R_trth[0] -6.95305e-13*Rs2_lac[6] +2.0e-10);


RF_s2_meanCorr_pad_7 = (Rs2_tdchit[15] - Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[23])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[55])*tdcTime)/2.0 + 0.802636e-08*R_trx[0]  -3.10881e-08*R_trth[0]-5.37148e-13*Rs2_lac[7]+ 6.200e-10 ); 
      
 RF_s2_meanCorr_pad_8 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[24])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[56])*tdcTime)/2.0 +0.787479e-08*R_trx[0]  - 3.18107e-08*R_trth[0]-5.95287e-13*Rs2_lac[8]); 

 RF_s2_meanCorr_pad_9 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[25])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[57])*tdcTime)/2.0 + 7.54862e-09*R_trx[0]  -3.5057e-08*R_trth[0] - 1.01789e-12*Rs2_lac[9] + 9.50e-10 ); 
      
RF_s2_meanCorr_pad_10 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[26])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[58])*tdcTime)/2.0 + 7.5127e-9*R_trx[0]  - 3.59703e-08*R_trth[0] - 1.02612e-12*Rs2_lac[10] + 1.0e-10);
    
RF_s2_meanCorr_pad_11 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[27])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[59])*tdcTime)/2.0 + 8.48865e-9*R_trx[0]  -3.76206e-08*R_trth[0] - 1.87664e-12*Rs2_lac[11]+ 2.2e-10 ); 
      
RF_s2_meanCorr_pad_12 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[28])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[60])*tdcTime)/2.0 + 5.42156e-9*R_trx[0]  - 3.68166e-08*R_trth[0] - 3.19282e-12*Rs2_lac[12] + 2.20e-9);


RF_s2_meanCorr_pad_13 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[29])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[61])*tdcTime)/2.0 + 6.27864e-9*R_trx[0]  - 3.51979e-08*R_trth[0] - 3.37812e-12*Rs2_lac[13]+ 2.0e-09 ); 


RF_s2_meanCorr_pad_14 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[30])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[62])*tdcTime)/2.0 + 5.78027e-9*R_trx[0]  - 3.54868e-08*R_trth[0] - 7.80362e-13*Rs2_lac[14] + 1.6e-09);

  RF_s2_meanCorr_pad_15 = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[31])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[63])*tdcTime)/2.0 + 8.605971e-10*R_trx[0]  - 1.56998e-08*R_trth[0] - 7.65658e-13*Rs2_lac[15]+ 1.615e-9 ); 
      
  double a[3] = {1,2,3};
  double b[3];
  b[0] = 1;
  b[1] =2;
  b[2] = 3;
  //a equals b
 

/////////////////////////////////////////////////////////////
 
  for( int i=2;i<16;i++)
    {
      if( HRS_trig_bit5 && Rs2_nthit==1 && Rs2_pads[0]==i )
	{
	  //Calc RF_S2mean depending on i
	  double RF_s2_mean = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime  - (((Rs2_tdchit[15]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0 + corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2]);
	    

// #some calculation depending on TDC information and correction depending on i ([paddle which was hit)
	  h1->Fill(RF_s2_mean);
	}
      else if ( HRS_trig_bit5 && Rs2_nthit==2 && Rs2_pads[0]==i && Rs2_pads[1]==i+1) 
	{
	}
      else {}
    }


  h1->Fill(RF_s2_mean);



/////////////////////////////////////////////////////////












			    










    } // end entry loop ie irun=90854 loop

    }// end subrun
    }

  TCanvas *c1 = new TCanvas("c1", "c1",600,600);
  c1->cd();
  h1->Draw();
}// end of void loop




