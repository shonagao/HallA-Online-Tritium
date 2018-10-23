//////////////////////////////////////////////////////////
//07/03/2018*!!!!!!!!!!!!!!!!!! Final code for the RHRS system. s2 detector
//Author Bishnu Pandey
//Here I am going to combine all of the scintillator together except 0 & 1
//First two scintillator doesnot have enough statistics, so they are not included in this code
///////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
///  FOr Marathan, runs start from 90854 to 90897 or so

/// For  coincidence, run starts from from 100400 to 100610


Int_t NEvent = 1000;

void RHRStime_1()
{
  gStyle->SetOptStat(1111);
   Double_t tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec
   TH2F *h1 = new TH2F("h1", "",500,0.30e-6,0.35e-6,17,1,16);
  
  TString filename;
  for(int irun= 100500;irun<100505;irun++)
    {

      if (irun == 90857) continue;
 

      for(int subrun =1; subrun<2; subrun++)
      	{
	 
      if( irun==90859 || subrun>>5) continue;

      //  filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d_%d.root", irun, subrun)); // for marathan run
       filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d.root", irun)); // For coincidence run
	   TFile *data_file = new TFile(filename,"READ");
           TTree *T = (TTree*)data_file->Get("T");
 
  // Variable to be used event by event 
	  Double_t RF_s2_mean;

   // Define Tree leaf variable to hold the values
            TString nHRS_trig = "DR.rrawt4"; // coincidence run
	  // TString nHRS_trig = "DR.evtypebits"; // for beta =1 particle
          TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
	   TString nRs2_tdchit = "RTDC.F1FirstHit"; // for coincidence run
	  // TString nRs2_tdchit = "F1FirstHit"; // for marathan run
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
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);


 Double_t corr_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};

	 
 Double_t corr_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};

	 
 Double_t corr_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};

	 
 Double_t alignment[14] = {2.0975e-9, 2.1e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
 
        /// loop over entries

	  Long64_t nentries = T->GetEntries();
	  for(Long64_t j=0;j<nentries;j++)
	    {
	      T->GetEntry(j); 
	 
     //UInt_t HRS_trig_bit5 =(static_cast<UInt_t>(HRS_trig)>>7)&1;  // For Coincidence run >>7)&1 and for Marathon >>5)&1
	 
              
/////////////////////////////////////////////////////////////
	      /// making loop for each paddle
  
for( int i=2;i<16;i++)
    {
      if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<0.80  ) // R.tr.beta>0.76 && R.tr.beta<80
	{
	  //Calc RF_S2mean depending on i
	    RF_s2_mean = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime  - (((Rs2_tdchit[15]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0 + corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2]);

	   
	      h1->Fill(RF_s2_mean, i);
	    }
// #some calculation depending on TDC information and correction depending on i ([paddle which was hit)
     	else if ( HRS_trig>1 && Rs2_nthit==2 && Rs2_pads[0]==i && Rs2_pads[1]==i+1 &&  R_trbeta>0.76 && R_trbeta<80) //R.tr.beta>0.76 && R.tr.beta<80
      	{
      	}
       else {}
    


    //  h1->Fill( RF_s2_mean,R_trth[0]);
    
////////////////////////////////////////////////////

// Here I am going to do same corrections for X and X' because this calibration is done for beta = 1 particle and in the RHRs system we do not have beta= 1 particle anymore.
  

  
 
////////////////////////////////////////////////////////
	    }



} // end entry loop ie irun=90854 loop

    }// end subrun
}

  TCanvas *c1 = new TCanvas("c1", "c1",600,600);
  c1->cd();
  h1->Draw();

 TLine *l2 = new TLine(0.525e-6,-0.2,0.525e-6,0.2);
 l2->SetLineColor(kRed);
  l2->Draw(); 


  //c1->SaveAs("RHRStime.pdf"); 

}// end of void loop

