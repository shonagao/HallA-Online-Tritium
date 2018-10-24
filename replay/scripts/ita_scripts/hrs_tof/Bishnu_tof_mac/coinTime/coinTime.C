
// Author Carlos and Bishnu
//06/07/2018/
// looking for an  event that either pass through pad#7 or pad# 8 or both



#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

void coinTime()
{
  
  Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch

  //Int_t runNUM = 90854;
  // Int_t evtNUM = 10000;


  
 TH1F *h7=new TH1F("h7","h7",500,0.55e-6,0.6e-6); 

      //TString filename = "/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_90854.root";
  //TFile *data_file = new TFile(filename, "READ");
  //TTree *T = (TTree*)data_file->Get("T");

  TString filename;

  for (int irun = 90854; irun < 90857; irun++)
    {
    
        if (irun==90858 ) continue;
	// if (irun==90860 ) continue;

      for (int subrun=1;subrun<9;subrun++)
	{

	  if((irun==90859 ||  subrun>>5) )continue;
	
	  // if(irun==90860 ||  subrun>>1)continue;
  
	  filename = (Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/Rootfiles/tritium_%d_%d.root", irun, subrun));
    
  TFile *data_file = new TFile(filename, "READ");
  TTree *T = (TTree*)data_file->Get("T");

  //  cout << "analyzing file: " << filename << endl; // This will find if any of the subrun root file is missing

 //Variables to be used event by event
  Double_t RF_s2_mean;
  Double_t RF_s2_pad7_meanCorr;
  Double_t RF_s2_pad8_meanCorr;


  //Define TTree Leaf names that will be used (shoule be exactly as it appears on the TTree)
  
  TString nHRS_trg = "DR.evtypebits";  //Ask how can you read DR.evtypebits>>5&1 to set as BranchAddress ???????
  TString nRs2_pad = "R.s2.t_pads";
  TString nRs2_nthit = "R.s2.nthit";
  TString nRs2_tdchit = "F1FirstHit";  
  TString nR_trx = "R.tr.x";
  TString nR_trth = "R.tr.th";
  TString nRs2_lac = "R.s2.la_c";

 

  //Define TTree Leaf variables to hold the values
  Double_t HRS_trg;
  Double_t Rs2_pad[100];  ///
  Double_t Rs2_nthit;
  Double_t Rs2_tdchit[100];  //need to find out the array size, use 100 for now  ??????
  //  Double_t* Rs2_tdchit = 0;
  Double_t  R_trx[100];
  Double_t R_trth[100];
  Double_t Rs2_lac[100];

  //Set the beanch address
  T->SetBranchAddress(nHRS_trg, &HRS_trg);
  T->SetBranchAddress(nRs2_pad, &Rs2_pad);
  T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
  T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
  T->SetBranchAddress(nR_trx, &R_trx);
  T->SetBranchAddress(nR_trth, &R_trth);
  T->SetBranchAddress(nRs2_lac,&Rs2_lac);

  //Loop over entries
  Long64_t nentries = T->GetEntries();
   
  for(Long64_t i=0; i<nentries; i++)
    {
      T->GetEntry(i);  

      UInt_t HRS_trg_bit5 = (static_cast<UInt_t>(HRS_trg)>>5)&1;

      //if (Rs2_pad[0]==7 && Rs2_pad[1]==8)
      //{
      
      //cout << "Event: " << i << endl;
      // cout << "HRS_trg = " << HRS_trg << endl;
      // cout << "HRS_trg_bit5 = " << HRS_trg_bit5 << endl;
      //cout << "Rs2_pad[0] = " << Rs2_pad[0] << endl;
      // cout << "Rs2_pad[1] = " << Rs2_pad[1] << endl;
      // cout << "Rs2nhit = " << Rs2_nthit << endl;
      
      /*
   for (int ch=0; ch < 64; ch++)
      	{
       cout << Form("Rs2_tdchit[%d] = ", ch) << Rs2_tdchit[ch] << endl;
	}      

       */
    
       RF_s2_pad7_meanCorr = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[23])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[55])*tdcTime)/2.0 + 0.802636e-08*R_trx[0]  -3.10881e-08*R_trth[0]-5.37148e-13*Rs2_lac[7]+ 6.200e-10 ); 
      

        // + 6.2e-10 is the best for #7
       // + 1.16319e-9  
       //  cout<< "RF Time ="<< RF_s2_pad7_meanCorr<<endl;
       //	h1->Fill(RF_s2_pad7_meanCorr);

 RF_s2_pad8_meanCorr = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime -(((Rs2_tdchit[15]-Rs2_tdchit[24])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[56])*tdcTime)/2.0 +0.787479e-08*R_trx[0]  - 3.18107e-08*R_trth[0]-5.95287e-13*Rs2_lac[8]);   
 //	h8->Fill(RF_s2_pad8_meanCorr);

	//Defining cuts 
	if( HRS_trg_bit5 && Rs2_pad[0]==7 &&  Rs2_nthit==1)
	  {
	    RF_s2_mean = RF_s2_pad7_meanCorr;
	  }
	else if( HRS_trg_bit5 && Rs2_pad[0]==8 && Rs2_nthit==1 )
	  {
               RF_s2_mean = RF_s2_pad8_meanCorr;
	  }

	else if(  HRS_trg_bit5 &&((Rs2_pad[0]==7 && Rs2_pad[1]==8) /*||(Rs2_pad[0]==7 || Rs2_pad[0]==8)*/))
	{
          RF_s2_mean = (RF_s2_pad7_meanCorr + RF_s2_pad8_meanCorr)/2.0;
	}
	else 
	{
	  RF_s2_mean = 0;
	}

	//hitso->Fille(RF_s2_mean);
        if (RF_s2_mean !=0)
	  {
	  h7->Fill(RF_s2_mean);
	  } 
	//	h8->Fill(RF_s2_mean);

    }  //end entry loop

	}//end subrun loop
 
    } //end run loop


  TCanvas *c1=new TCanvas("c1","c1",1200,1200);
  c1->cd();
  h7->Draw();
  // TCanvas *c8 = new TCanvas("c8","c8",600,600);
  //// c8->cd();
  // h8->Draw();
}



 
