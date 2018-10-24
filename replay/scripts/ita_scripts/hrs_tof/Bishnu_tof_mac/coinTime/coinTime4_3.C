// Author  Bishnu
//06/13/2018/
// looking for an  event that either pass through pad#4 or pad# 4 or both
// coping from coinTime6_5.C

#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

void coinTime4_3()
{
 Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch

 TH1F *h3=new TH1F("h3","either#3 hit or #4 or both;RF -s2 time in sec;counts",500,0.55e-6,0.6e-6);
 TString filename;
 for(int irun=90854;irun<90858;irun++)
   {
     if(irun==90858) continue;
     for(int subrun=1;subrun<9;subrun++)
       {
	 if(irun==90859 ||subrun>>5) continue;

     filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d_%d.root", irun, subrun));

     TFile *data_file = new TFile(filename,"READ");
     TTree *T = (TTree*)data_file->Get("T");

     //variable to be used  event by event
  Double_t RF_s2_mean;
 Double_t RF_s2_pad3_meanCorr;
  Double_t RF_s2_pad4_meanCorr;
  
  // Define TTree leaf name that will be used(should be exactly same as it appear on the TTree)
  TString nHRS_trig = "DR.evtypebits";
  TString nRs2_pads = "R.s2.t_pads";
  TString  nRs2_nthit = "R.s2.nthit";
  TString nRs2_tdchit = "F1FirstHit";
  TString nR_trx = "R.tr.x";
  TString nR_trth = "R.tr.th";
  TString nRs2_lac = "R.s2.la_c";


//// Define Tree leaf variables to hold the values
  Double_t HRS_trig;
  Double_t Rs2_pads[100];
  Double_t Rs2_nthit;
  Double_t Rs2_tdchit[100]; // need to find the aray size 100 for now
  Double_t R_trx[100];
  Double_t R_trth[100];
  Double_t Rs2_lac[100];

  // Set the branch address
  T->SetBranchAddress(nHRS_trig, &HRS_trig);
  T->SetBranchAddress(nRs2_pads, &Rs2_pads);
  T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
  T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
  T->SetBranchAddress(nR_trx,&R_trx);
  T->SetBranchAddress(nR_trth,&R_trth);
  T->SetBranchAddress(nRs2_lac, &Rs2_lac);




  ///Loop over entries
  Long64_t nentries = T->GetEntries();
  for(Long64_t i=0;i<nentries;i++)
    {
      T->GetEntry(i);
      UInt_t HRS_trig_bit5 = (static_cast<UInt_t>(HRS_trig)>>5)&1;
      if(HRS_trig_bit5 ==0) continue;

RF_s2_pad3_meanCorr = (Rs2_tdchit[15] - Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[19])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[51])*tdcTime)/2.0 +5.41034e-09 *R_trx[0] -4.79507e-08*R_trth[0] - 3.84824e-13*Rs2_lac[5] + 2.1e-9); 
// 2.1 e-9 gives , 336,353 and 360 ps


RF_s2_pad4_meanCorr = (Rs2_tdchit[15]-Rs2_tdchit[14])*tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[20])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[52])*tdcTime)/2.0 + 6.87688e-09*R_trx[0] -3.59540e-08*R_trth[0] -1.45016e-12*Rs2_lac[4] +0.85e-9 );



 


			      // Defining cuts


			      if(HRS_trig_bit5 && Rs2_pads[0]==3 && Rs2_nthit==1)
				{
				  RF_s2_mean =  RF_s2_pad3_meanCorr;
				}
			     else if( HRS_trig_bit5 && Rs2_pads[0]==4 && Rs2_nthit==1 )
				{
                                RF_s2_mean =  RF_s2_pad4_meanCorr;
				}
			     else if( HRS_trig_bit5 && Rs2_pads[0]==3 && Rs2_pads[1]==4 )
			       {

                                 RF_s2_mean =(RF_s2_pad3_meanCorr + RF_s2_pad4_meanCorr)/2.0;

			       }
			     else
			       {
				 RF_s2_mean =0;
			       }



			      if(RF_s2_mean!=0)
				{
				  h3->Fill(RF_s2_mean);
				}
    }// end entry loop






     // cout << "dta file= "<< data_file<<endl;


     // cout<<"File running ="<<filename<<endl;



       } // end subrun loop
   } // end irun  loop

 TCanvas *c3=new TCanvas("c3","c3",1200,1200);
  c3->cd();
  h3->Draw();
}// void
