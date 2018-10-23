/// combining paddle #13 and #13

#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

void L_add13_14(){
  
  Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch

  //Int_t runNUM = 90854;
  // Int_t evtNUM = 10000;
TH1F *h14=new TH1F("h14","h14 with  +7.9375e-10 ",500,0.10e-6,0.15e-6);
TString filename;
 for(int irun=100380;irun<100585;irun++) {

   if(irun!=100397 && irun!=100398 && irun!=100399 && irun!= 100422 && irun!=100423 && irun!=100441 && irun!=100442 && irun!= 100446 && irun !=100450 && irun!=100451 && irun!=100460 && irun!=100477 && irun!=100478 && irun!=100479 && irun!=100480 && irun!=100492 && irun!=100510 && irun!=100511 && irun!=100522 && irun!=100523 && irun!=100524 && irun!=100541 && irun!=100549 && irun!=100550 && irun!=100555 && irun!=100557 && irun!=100563 && irun!=100567 )
     {
     for (int subrun = 1;subrun<2;subrun++)
       {

	 filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d_%d.root",irun,subrun));
	 TFile *data_file = new TFile(filename,"READ"); 
	 TTree *T = (TTree*)data_file->Get("T");

	 /// variables to be used event by event

     Double_t RF_s2_mean;
     Double_t RF_s2_pad13_meanCorr;
     Double_t RF_s2_pad14_meanCorr;

//Define TTree Leaf names that will be used (shoule be exactly as it appears on the TTree)
     TString nLHRS_trig = "DR.evtypebits";
     TString nLs2_pad = "L.s2.t_pads";
     TString nLs2_nthit = "L.s2.nthit";
     TString nLs2_tdchit = "LTDC.F1FirstHit";
     TString nL_trx = "L.tr.x";
     TString nL_trth = "L.tr.th";
     TString nLs2_ladc = "L.s2.la_c";

// Define TTree Leaf  variables to hold the values
    Double_t LHRS_trig;
    Double_t  Ls2_pad[100];
    Double_t Ls2_nthit;
    Double_t Ls2_tdchit[100];
    Double_t L_trx[100];
    Double_t  L_trth[100];
    Double_t Ls2_ladc[100];

 // Now set the branch address
    T->SetBranchAddress(nLHRS_trig, &LHRS_trig);
    T->SetBranchAddress(nLs2_pad, &Ls2_pad);
    T->SetBranchAddress(nLs2_nthit, &Ls2_nthit);
    T->SetBranchAddress(nLs2_tdchit, &Ls2_tdchit);
    T->SetBranchAddress(nL_trx, &L_trx);
    T->SetBranchAddress(nL_trth, &L_trth);
    T->SetBranchAddress(nLs2_ladc, &Ls2_ladc);


		     // loop over nentries
		     Long64_t nentries = T->GetEntries();
		     for(Long64_t i=0;i<nentries;i++)
		       {
			 T->GetEntry(i);
			 UInt_t LHRS_trig_bit2 = (static_cast<UInt_t>(LHRS_trig)>>2)&1;

			 // Do not use these variables
 
   





 RF_s2_pad13_meanCorr =(Ls2_tdchit[37] -Ls2_tdchit[47])*tdcTime  - (((Ls2_tdchit[30]- Ls2_tdchit[13])*tdcTime + (Ls2_tdchit[37] - Ls2_tdchit[61])*tdcTime )/2.0 + 6.22836e-09*L_trx[0]  -3.07010e-08*L_trth[0] -1.23475e-12*Ls2_ladc[13] + 1.485e-10); 

RF_s2_pad14_meanCorr = (Ls2_tdchit[37] -Ls2_tdchit[47])*tdcTime  - (((Ls2_tdchit[30] - Ls2_tdchit[14])*tdcTime + (Ls2_tdchit[37] - Ls2_tdchit[62])*tdcTime )/2.0 + 5.52765e-09*L_trx[0] -3.79624e-08*L_trth[0] -1.50406e-12*Ls2_ladc[14] + 7.9375e-10);

 


    

 // Defining cuts
	   if(LHRS_trig_bit2 && Ls2_pad[0]==13 && Ls2_nthit==1)
	   {
             RF_s2_mean = RF_s2_pad13_meanCorr;
	   }
	   else if( LHRS_trig_bit2 && Ls2_pad[0]==14 && Ls2_nthit==1 )
	     {
                RF_s2_mean = RF_s2_pad14_meanCorr;
	     }
	   else if( LHRS_trig_bit2 && (Ls2_pad[0]==13 && Ls2_pad[1]==14 ))

	     {
            RF_s2_mean = (RF_s2_pad13_meanCorr + RF_s2_pad14_meanCorr)/2.0;
	}

	   else
		 {
                  RF_s2_mean = 0;
		 }

                if (RF_s2_mean !=0)
	  {
	    h14->Fill(RF_s2_mean);
	  } 



	 } // end Entry loop

       }// end of subrun
   } // end run loop
 }
 TCanvas *c1 = new TCanvas("c1","c1",600,600);
 c1->cd();
 h14->Draw();


}
