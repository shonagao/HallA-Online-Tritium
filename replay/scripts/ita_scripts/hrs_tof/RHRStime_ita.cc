//////////////////////////////////////////////////////////
//07/03/2018*!!!!!!!!!!!!!!!!!! Final code for the RHRS system. s2 detector
//Author Bishnu Pandey
//Here I am going to combine all of the scintillator together except 0 & 1
//First two scintillator doesnot have enough statistics, so they are not included in this code
///////////////////////////////////////////////////////////////////
// Here I am corecting to the X and X' for proton by using the coincidence run
//here I used the coincidence runs or detector is now corrected for(ee'p) runs
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
///  FOr Marathan, runs start from 90854 to 90897 or so

/// For  coincidence, run starts from from 100400 to 100610
// correcting for X and X' again for the ee'p data//////07/24/2018
// just copied from RHRStime_3.C as RHRStime_4.C doesn't look good
// just for the testing purpose. probably have to do the X, X' corrections .. 10/01/2018
// This code I am using R1A = 8 and RF time = 1 for few days
Int_t NEvent = 1000;

void RHRStime_ita()
{
  gStyle->SetOptStat(111111);
   Double_t tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec
   // TH2F *h1 = new TH2F("h1", "",700,0.30e-6,0.35e-6,300,-50,50); //For Two Dimensional Histogram
   //  TH2F *h1 = new TH2F("h1", "",500,0.3e-6,0.35e-6,14,0.5,15.5);// one dimensional Histogram 
  TH1F *h1 = new TH1F("h1", "",1000,0.65e-6,0.7e-6);// one dimensional Histogram
   //  TProfile *p1 = new TProfile("p1", "",200,-1,1,0.3190e-6,0.3221e-6);
  TString filename;
  for(int irun= 93498;irun<93510;irun++)
    {
      if (irun == 100422) continue;
  for(int subrun =1; subrun<2; subrun++)
      	{	 
      if( irun==100422 || subrun>>5) continue;

      // filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d_%d.root", irun, subrun)); // for marathan run
      filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d.root", irun)); // For coincidence run
	   TFile *data_file = new TFile(filename,"READ");
           TTree *T = (TTree*)data_file->Get("T");
 
  // Variable to be used event by event 
	  Double_t RF_s2_mean;
   // Define Tree leaf variable to hold the values
            TString nHRS_trig = "DR.rrawt4";
	  //TString nHRS_trig = "DR.evtypebits"; // for beta =1 particle
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

	  Long64_t nentries = T->GetEntries();
	  for(Long64_t j=0;j<nentries;j++)
	    {
	      T->GetEntry(j); 
	     
/////////////////////////////////////////////////////////////
	      /// making loop for each paddle
  
for( int i=2;i<16;i++)
    {
      if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01  ) // R.tr.beta>0.76 && R.tr.beta<80
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

}// end of void loop

