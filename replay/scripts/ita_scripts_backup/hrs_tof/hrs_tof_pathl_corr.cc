#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"


Int_t NEvent = 1000;

void hrs_tof_pathl_corr()
{
  gStyle->SetOptStat(1111);
    Double_t tdcTime = 56.23e-12; // this is the conversion per chanel ie 56.23 picosec
    TH2F *h1 = new TH2F("h1", "",300,-0.11,0.11,500,0.5e-6,0.52e-6); //For Two Dimensional Histogram
    TH2F *h2 = new TH2F("h2", "",300,-1,1,500,0.5e-6,0.52e-6); //For Two Dimensional Histogram
    // TH2F *htof_pathl = new TH2F("htof_pathl", "",100,0.5e-6,0.52e-6,50,22,23); //For Two Dimensional Histogram

 TH2F *htof_pathl = new TH2F("htof_pathl", "",10000,-1.0e-5,1.0e-6,50,22,23); //For Two Dimensional Histogram
   //  TH2F *h1 = new TH2F("h1", "",500,0.3e-6,0.35e-6,14,0.5,15.5);// one dimensional Histogram 
   // TH1F *h1 = new TH1F("h1", "",1000,0.65e-6,0.7e-6);// one dimensional Histogram
    TProfile *p1 = new TProfile("p1", "",100,-0.11,0.11,0.5093e-6,0.5113e-6); //For Two Dimensional Histogram
    TProfile *p2 = new TProfile("h2", "",100,-1,1,0.5093e-6,0.5113e-6); //For
    // Two Dimensional Histogram
 



 TChain*  T=new TChain("T");
  for(int irun= 94010;irun<94017;irun++) // 93498 tp 93510 for the short range run // 93774 October 8 run
    {
T->Add (Form("/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root", irun)); // For coincidence run 
    }

	   
  // Variable to be used event by event 
	  Double_t RF_s2_mean;

   // Define Tree leaf variable to hold the values
	  TString nHRS_trig = "DR.T5";
	  // TString nHRS_trig = "DR.rrawt4";
	  //TString nHRS_trig = "DR.evtypebits"; // for beta =1 particle
          TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
	  TString nRs2_tdchit = "RTDC.F1FirstHit"; // for coincidence run
	  // TString nRs2_tdchit = "F1FirstHit"; // for marathan run
	  TString nR_trx = "R.tr.x";
          TString nR_trth = "R.tr.th";
          TString nRs2_lac = "R.s2.la_c";
          TString nR_trbeta="R.tr.beta";
	   TString Rpathl_s2="R.s2.trpath";
	   TString Rpathl_rf="R.tr.pathl";

    // Define Tree leaf variable
	  Double_t HRS_trig;
	  Double_t Rs2_pads[100];
	  Double_t Rs2_nthit;
	  Double_t Rs2_tdchit[100];
	  Double_t R_trx[100];
	  Double_t R_trth[100];
	  Double_t Rs2_lac[100];
          Double_t R_trbeta;
	  Double_t pathl_s2[1000];
	  Double_t pathl_rf[1000];
    // Set The Branch addres
          T->SetBranchStatus("*",0);
	  T->SetBranchStatus(nHRS_trig,1);
          T->SetBranchStatus(nRs2_pads,1);
          T->SetBranchStatus(nRs2_nthit,1);
          T->SetBranchStatus(nRs2_tdchit,1);
          T->SetBranchStatus(nR_trx,1);
          T->SetBranchStatus(nR_trth,1);
          T->SetBranchStatus(nRs2_lac,1);
          T->SetBranchStatus(nR_trbeta,1);
	  T->SetBranchStatus(Rpathl_s2,1);
	  T->SetBranchStatus(Rpathl_rf,1);
         
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);
	  T->SetBranchAddress(Rpathl_s2,pathl_s2);
	  T->SetBranchAddress(Rpathl_rf,pathl_rf);

	  Double_t pathl;	
    // original values
 Double_t corr_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};
 /// the following loop changes from Marathon to coincidence runs
 
// Double_t c =  -3.9094e-09;  // -3.01594e-9 // real value -3.9094e-09; 

 Double_t c_x =  -1.584e-10;  // -3.01594e-9 // real value -3.9094e-09; 
 for(int i=0;i<14;i++)
   {
     corr_x[i] = corr_x[i] + c_x;
   }

 
 // since I am using the coincidence data now and there is some offset. I want to remove the offset. So I would like to add a constant in each of the array element.	 
 Double_t corr_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};
 
 Double_t c_th =  -2.075e-09;  // -3.01594e-9 // real value -3.9094e-09; 
 for(int i=0;i<14;i++)
   {
     corr_th[i] = corr_th[i] + c_th;
   }
 
 

	 
 Double_t corr_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};
 // the following corrections will not be applied

 	 
 // Double_t alignment[14] = {2.0975e-9, 2.1e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
   Double_t alignment[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
   // This value is taken from coin_trig.C which is the common trig for LHRS and RHrs. Actually I took this value from test_coin.C
        /// loop over entries

int ch=8;

	  Long64_t nentries = T->GetEntries();
	  cout<<"Get Entries : "<<nentries<<endl;
	  for(Long64_t j=0;j<nentries;j++){
	      T->GetEntry(j);    
	    for( int i=2;i<16;i++){
              if(i==ch){
       
		//         if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==i &&  R_trbeta>0.76 && R_trbeta<1.01  ){ // R.tr.beta>0.76 && R.tr.beta<80	
               
	   
	     RF_s2_mean = (Rs2_tdchit[9] - Rs2_tdchit[15]) * tdcTime  - (((Rs2_tdchit[9]-Rs2_tdchit[i+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[i+48])*tdcTime)/2.0 + corr_x[i-2]*R_trx[0] + corr_th[i-2]*R_trth[0] + corr_adc[i-2]*Rs2_lac[i] + alignment[i-2]); // late ST data
	    
	     pathl =pathl_rf[0]-pathl_s2[0];
	   
	  //  h1->Fill(RF_s2_mean);// for 2 Dimensional;

	  h1->Fill(R_trth[0], RF_s2_mean);// For opne  dimensional
	  p1->Fill(R_trth[0], RF_s2_mean);// For opne  dimensional
	  h2->Fill(R_trx[0], RF_s2_mean);// For opne  dimensional
	  p2->Fill(R_trx[0], RF_s2_mean);// For opne  dimensional
	  htof_pathl->Fill(RF_s2_mean,pathl);
	 
	  // }
      
  }
	      }

	  } // end entry loop ie irun=90854 loop


  TF1 *f1 = new TF1("f1","pol1",-0.11,0.11);
  p1->Fit(f1,"0","",-0.11,0.11);

  // TF1 *f1 = new iTF1("f1","[0]+[1]*x",-5,12);
  // h1->Fit("f1","R+");
 
  TCanvas *c1 = new TCanvas("c1", "c1",600,600);
  c1->cd();
  h1->Draw("colz");
  p1->Draw("same");
  f1->Draw("same");
  TCanvas *c2 = new TCanvas("c2", "c2",600,600);
  c2->cd();
  h2->Draw("colz");
  p2->Draw("same");
  TCanvas *c3 = new TCanvas("c3", "c3",600,600);
  c3->cd();
  h2->ProjectionY()->Draw("");

  TCanvas *ctof_pathl = new TCanvas("ctof_pathl", "ctfo_pathl",600,600);
  ctof_pathl->cd();
  htof_pathl->Draw("colz");

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

