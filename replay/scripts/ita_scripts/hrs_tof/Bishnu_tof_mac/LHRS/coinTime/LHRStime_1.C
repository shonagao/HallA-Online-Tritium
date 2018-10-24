//////////////////////////////////////////////////////////
//07/19/2018******** Final code for the LHRS system.
//Author Bishnu Pandey
//Here I am going to combine all of the scintillator together except 0 & 15
//First and last scintillator doesnot have enough statistics, so they are not included in this code
///////////////////////////////////////////////////////////////////

///Note the pmt are not aligned well. need some more work.
// making Rf -s2 vs L.trx and L.s2.t_pads


#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
void LHRStime_1()
{
  
  Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch

  //Int_t runNUM = 90854;
  // Int_t evtNUM = 10000;
  TH2F *h1=new TH2F("h1"," LHRS ",500,0.10e-6,0.155e-6,14,0,16);
 TString filename;
 for(int irun=100400;irun<100401;irun++) //100380;irun<100585
   {
 if(irun!=100397 && irun!=100398 && irun!=100399 && irun!= 100422 && irun!=100423 && irun!=100441 && irun!=100442 && irun!= 100446 && irun !=100450 && irun!=100451 && irun!=100460 && irun!=100477 && irun!=100478 && irun!=100479 && irun!=100480 && irun!=100492 && irun!=100510 && irun!=100511 && irun!=100522 && irun!=100523 && irun!=100524 && irun!=100541 && irun!=100549 && irun!=100550 && irun!=100555 && irun!=100557 && irun!=100563 && irun!=100567 )
     {

     for (int subrun = 1;subrun<2;subrun++)

       {

	 filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d_%d.root",irun,subrun));
	 TFile *data_file = new TFile(filename,"READ"); 
	 TTree *T = (TTree*)data_file->Get("T");

	 /// variables to be used event by event
  Double_t RF_s2_mean;
 
 


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

 Double_t corr_x[14] = {9.5982e-09, 2.39686e-09, 5.50452e-09, 8.67284e-09, 7.88134e-09, 9.39930e-09,   9.09441e-09, 8.13305e-09, 8.36477e-09, 8.74297e-09, 7.745e-09,  5.94972e-09, 6.22836e-09, 5.52765e-09};

 Double_t corr_th[14] = {-5.3783e-08,  - 3.32898e-08, -4.14532e-08, -4.08767e-08, -4.07972e-08, -3.63437e-08,  -3.67840e-08, -3.54952e-08, -3.63706e-08,-3.39145e-08, -3.43925e-08,  -3.05521e-08,-3.07010e-08, -3.79624e-08};


Double_t corr_adc[14] = {- 1.592e-12, -1.24122e-12, -1.18518e-12, -1.16133e-12, -1.24632e-12, -1.22617e-12, -1.02470e-12, -6.57058e-13, -1.14584e-12, -1.3259e-12, -1.816135e-12, -1.15547e-12,  -1.23475e-12, -1.50406e-12};

// Double_t alignment[14] = {3.19760e-10, 3.078e-9, 3.435e-9, 9.985e-10, 9.835e-10,  4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10               };
 Double_t alignment[14] = {1.0319760e-9, -1.0e-9, -0.35e-9, 9.985e-10, 9.835e-10,  4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10 }; /// This came fromcoin_trig.C which is the coincidence trigger for LHRS and RHRS. Actually thei scame from test_coin.C



/// loop over entries
 Long64_t nentries = T->GetEntries();
		     for(Long64_t i=0;i<nentries;i++)
		       {
			 T->GetEntry(i);

			 UInt_t LHRS_trig_bit2 = (static_cast<UInt_t>(LHRS_trig)>>2)&1;

			 if(LHRS_trig_bit2 ==0)  continue;
			 //making loop fo reach paddle

			 for (int i=1;i<15;i++)
			   {
			     if(LHRS_trig_bit2 && Ls2_nthit==1 &&  Ls2_pad[0]==i)
                            {
			       RF_s2_mean = (Ls2_tdchit[37] -Ls2_tdchit[47])*tdcTime - (((Ls2_tdchit[30] - Ls2_tdchit[i])*tdcTime + (Ls2_tdchit[37] - Ls2_tdchit[i+48])*tdcTime )/2.0 + corr_x[i-1]*L_trx[0] + corr_th[i-1]*L_trth[0] + corr_adc[i-1]*Ls2_ladc[i] + alignment[i-1]);



			    
			      h1->Fill(RF_s2_mean,i);
			    }
			     else if(LHRS_trig_bit2 && Ls2_nthit==2&& Ls2_pad[0]==i && Ls2_pad[1]==i+1)
			       {
			       }
			     else{}
			     h1->Fill(RF_s2_mean,i);
			   }		
		       } // end of entry loop
       } /// end of sub run

     }
   }
TCanvas *c1 = new TCanvas("c1","c1",600,600);
 c1->cd();

 h1->Draw("colz");
 TLine *l1 = new TLine(0.125e-6,0,0.125e-6,16);
 l1->SetLineColor(kRed);
  l1->Draw(); 

   
}
