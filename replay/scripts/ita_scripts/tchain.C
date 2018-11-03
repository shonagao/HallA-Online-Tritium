//November 02,2018
// Looking for the ee'K data
//////////////////////////////////////////////////////////////////
/// Author Bishnu Pandey   07/26/2018
////This is the coincidence time for Left and right arms.
//////////////////////////////////////////////////////////////////////
// This code is working for the early ee'k data.
//here I am using TChain method to include all of root files ie subrootfiles
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
void tchain()
{

 TChain *T = new TChain("T");
  Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch
  gStyle->SetOptStat(11111);
  //  TH2F *h1 = new TH2F("h1","", 500,-0.145e-6,-0.12e-6,14,1,15);
  TH1F *h1 = new TH1F("h1","", 4000,-15,15); // -0.145e-6,-0.12e-6);
  //  TString filename;
  for(int irun = 111190;irun<111200;irun++)//100677;irun<100679
    {
      if(irun ==111193)continue;
      T->Add("/adaqfs/home/a-onl/tritium/replay/t2root/Rootfiles/tritium_111170*.root");
	 
	  //$$ filename = (Form("/adaqfs/home/a-onl/tritium/replay/t2root/Rootfiles/tritium_%d.root", irun));
      // TFile *data_file = new TFile(filename,"Read");
      //TTree *T = (TTree*)data_file->Get("T");

	  Double_t s2L;
	  Double_t s2R;
	  Double_t coin_trig;
	  /////////////////////////////////////////////////////////
	  ///Defining variables for Left HRS 
	  /////////////////////////////////////////////////////////

	  /// variables to be used for Left HRAS
	  Double_t RF_s2L_mean;
	  // Define TTree leaf names for Left Hrs that will be used.
	  // TString nLHRS_trig = "DR.evtypebits";
	  TString nLHRS_trig = "DR.T1";
	  TString nLs2_pad = "L.s2.t_pads";
	  TString nLs2_nthit = "L.s2.nthit";
	  TString nLs2_tdchit = "LTDC.F1FirstHit";
	  TString nL_trx = "L.tr.x";
	  TString nL_trth = "L.tr.th";
	  TString nLs2_ladc = "L.s2.la_c";
	  // Define TTree leaf variables to hold the values
	  Double_t LHRS_trig;
	  Double_t  Ls2_pad[100];
	  Double_t Ls2_nthit;
	  Double_t Ls2_tdchit[100];
	  Double_t L_trx[100];
	  Double_t  L_trth[100];
	  Double_t Ls2_ladc[100];

	  // Now set the branch address LHRS
	  T->SetBranchAddress(nLHRS_trig, &LHRS_trig);
	  T->SetBranchAddress(nLs2_pad, &Ls2_pad);
	  T->SetBranchAddress(nLs2_nthit, &Ls2_nthit);
	  T->SetBranchAddress(nLs2_tdchit, &Ls2_tdchit);
	  T->SetBranchAddress(nL_trx, &L_trx);
	  T->SetBranchAddress(nL_trth, &L_trth);
	  T->SetBranchAddress(nLs2_ladc, &Ls2_ladc);
	  // Defining array LHRS
	  Double_t corr_L_x[14] = {9.5982e-09, 2.39686e-09, 5.50452e-09, 8.67284e-09, 7.88134e-09, 9.39930e-09,   9.09441e-09, 8.13305e-09, 8.36477e-09, 8.74297e-09, 7.745e-09,  5.94972e-09, 6.22836e-09, 5.52765e-09};
	  
	  Double_t cL =-9.15734e-10;   // previous one  4.87486E-11;                -3.9094e-09;
	  for(int l=0;l<14;l++)
	    {
	      corr_L_x[l] = corr_L_x[l] + cL;
	    }
	  



	  Double_t corr_L_th[14] = {-5.3783e-08,  - 3.32898e-08, -4.14532e-08, -4.08767e-08, -4.07972e-08, -3.63437e-08,  -3.67840e-08, -3.54952e-08, -3.63706e-08,-3.39145e-08, -3.43925e-08,  -3.05521e-08,-3.07010e-08, -3.79624e-08};


	  Double_t cL1 = 1.75759e-9;   // previous one  4.87486E-11;                -3.9094e-09;
	  for(int m=0;m<14;m++)
	    {
	      corr_L_th[m] = corr_L_th[m] + cL1;
	    }
	  




	  
	  Double_t corr_L_adc[14] = {- 1.592e-12, -1.24122e-12, -1.18518e-12, -1.16133e-12, -1.24632e-12, -1.22617e-12, -1.02470e-12, -6.57058e-13, -1.14584e-12, -1.3259e-12, -1.816135e-12, -1.15547e-12,  -1.23475e-12, -1.50406e-12};

 Double_t alignment_L[14] = {1.0319760e-9, -1.0e-9, -0.35e-9, 9.985e-10, 9.835e-10,  4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10 };

 //Double_t alignment_L[14] = {7.319760e-9, -1.0e-9, -0.35e-9, 9.985e-10, 9.835e-10,  4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10 };


	
	  /////////////////////////////////////////////
	  //now defining same for the RHRS system
	  ////////////////////////////////////////


	  // Defining variables for RHRS 

	  Double_t RF_s2R_mean;
	  // Define Tree Leaf name for the RHRS that will be used.

	  //TString nHRS_trig = "DR.rrawt4";
	  TString nHRS_trig = "DR.T4";
	  //TString nHRS_trig = "DR.evtypebits";
	  TString nRs2_pads = "R.s2.t_pads";
          TString nRs2_nthit = "R.s2.nthit";
	  TString nRs2_tdchit = "RTDC.F1FirstHit";
	  TString nR_trx = "R.tr.x";
          TString nR_trth = "R.tr.th";
          TString nRs2_lac = "R.s2.la_c";
          TString nR_trbeta="R.tr.beta";
	  TString nR_a1 ="R.a1.asum_c";
	  TString nR_a2 ="R.a2.asum_c";
	  // NOw define the tree leaf variable
	  Double_t HRS_trig;
	  Double_t Rs2_pads[100];
	  Double_t Rs2_nthit;
	  Double_t Rs2_tdchit[100];
	  Double_t R_trx[100];
	  Double_t R_trth[100];
	  Double_t Rs2_lac[100];
          Double_t R_trbeta;
	  Double_t R_a1;
	  Double_t R_a2;

	  // Set the branch address
	  T->SetBranchAddress(nHRS_trig, &HRS_trig);
          T->SetBranchAddress(nRs2_pads, &Rs2_pads);
          T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
          T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
          T->SetBranchAddress(nR_trx, &R_trx);
          T->SetBranchAddress(nR_trth, &R_trth);
          T->SetBranchAddress(nRs2_lac, &Rs2_lac);
          T->SetBranchAddress(nR_trbeta, &R_trbeta);
	  T->SetBranchAddress(nR_a1, &R_a1);
	  T->SetBranchAddress(nR_a2, &R_a2);


	  // Defining some array variables

	  Double_t corr_R_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};

          
	  Double_t c =4.87486E-11;   // previous one  4.87486E-11;                -3.9094e-09;
	  for(int l=0;l<14;l++)
	    {
	      corr_R_x[l] = corr_R_x[l] + c;
	    }
	  
	  

	  Double_t corr_R_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};
	  //-3.06204e-09; 
	  
	  Double_t c1 =-3.06204E-09 ;   // previous one-3.06204E-09;     -3.9094e-09;
	  for(int m=0;m<14;m++)
	    {
	      corr_R_th[m] = corr_R_th[m] + c1;
	    }
	  
 

	  Double_t corr_R_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13,-6.95305e-13, -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};


	  // Double_t alignment_R[14] = {2.0975e-9, 2.1e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
	  //Double_t alignment_R[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 2.0e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};

 Double_t alignment_R[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};


	  
	  // Loop over entries

	  Long64_t nentries = T->GetEntries();
	  //	  nentries = 10000;
	  for(Long64_t k=0;k<nentries;k++)
	    {
	      T->GetEntry(k); 
	      /////////////////////////////////////////////////
	      ///// loop for LHRS
	      /////////////////////////////////////////////////////////////////////////////////////
	      
	      UInt_t LHRS_trig_bit2 = (static_cast<UInt_t>(LHRS_trig)>>5)&1; // was 2 before

	      // if(LHRS_trig_bit2 ==0)  continue;

			 for (int i=1;i<15;i++)
			   {
			     if(LHRS_trig_bit2 && Ls2_nthit==1 &&  Ls2_pad[0]==i)
                            {
			      // RF_s2L_mean = (Ls2_tdchit[37] -Ls2_tdchit[47])*tdcTime - (((Ls2_tdchit[30] - Ls2_tdchit[i])*tdcTime + (Ls2_tdchit[37] - Ls2_tdchit[i+48])*tdcTime )/2.0 + corr_L_x[i-1]*L_trx[0] + corr_L_th[i-1]*L_trth[0] + corr_L_adc[i-1]*Ls2_ladc[i] + alignment_L[i-1]);	// original 	     
RF_s2L_mean = (((Ls2_tdchit[30] - Ls2_tdchit[i])*tdcTime + (Ls2_tdchit[37] - Ls2_tdchit[i+48])*tdcTime )/2.0 + corr_L_x[i-1]*L_trx[0] + corr_L_th[i-1]*L_trth[0] + corr_L_adc[i-1]*Ls2_ladc[i] + alignment_L[i-1]);	

    ///////////////////////////////////////////////////////////////////////////////////
			 //Loop for RHRS paddles
	      //////////////////////////////////////////////////////////////////////
	      
	      for( int j=2;j<16;j++)
		{
		  if( HRS_trig>1 && Rs2_nthit==1 && Rs2_pads[0]==j && R_trbeta>0.76 && R_trbeta<2 && R_a1<150.0 && R_a2>1400.0 ) // R.a1.asum_p<150. && R.a2.asum_p>1400
		    {
		      //Calc RF_S2mean depending on i
		      // RF_s2R_mean = (Rs2_tdchit[9] - Rs2_tdchit[15]) * tdcTime  - (((Rs2_tdchit[9]-Rs2_tdchit[j+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[j+48])*tdcTime)/2.0 + corr_R_x[j-2]*R_trx[0] + corr_R_th[j-2]*R_trth[0] + corr_R_adc[j-2]*Rs2_lac[j] + alignment_R[j-2]);	   // original
	RF_s2R_mean =  (((Rs2_tdchit[9]-Rs2_tdchit[j+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[j+48])*tdcTime)/2.0 + corr_R_x[j-2]*R_trx[0] + corr_R_th[j-2]*R_trth[0] + corr_R_adc[j-2]*Rs2_lac[j] + alignment_R[j-2]);		//  RF_s2R_mean = (Rs2_tdchit[15] - Rs2_tdchit[14]) * tdcTime  - (((Rs2_tdchit[15]-Rs2_tdchit[j+16])*tdcTime + (Rs2_tdchit[46]-Rs2_tdchit[j+48])*tdcTime)/2.0 + corr_R_x[j-2]*R_trx[0] + corr_R_th[j-2]*R_trth[0] + corr_R_adc[j-2]*Rs2_lac[j] + alignment_R[j-2]);	     
                       
		       coin_trig = (RF_s2L_mean -RF_s2R_mean);
		       coin_trig = coin_trig*1e+9 -245.832;      /// changing in tonano second scale
		       h1->Fill(coin_trig);
		       //    cout << coin_trig << endl;
		    
		     		    }
		    //  else if ( HRS_trig>1 && Rs2_nthit==2 && Rs2_pads[0]==i && Rs2_pads[1]==j+1 &&  R_trbeta>0.76 && R_trbeta<80) //  R.tr.beta>0.76 && R.tr.beta<80
			  //   {
		    //  }
		    // else {}		 
		} // end of j<16 loop
                              }			     
	      	      else if(LHRS_trig_bit2 && Ls2_nthit==2&& Ls2_pad[0]==i && Ls2_pad[1]==i+1)
			       {
			        }
	     	      else{}
	       }/// end of i<15 loop
		
  

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////	  

	   }// end of entry loop


	

     
    }// end of irun loop

 TCanvas *c1 = new TCanvas("c1","c1",600,600);
 c1->cd();
 h1->Draw("colz");
  TLine *l1 = new TLine(-0.132e-6,1,-0.132e-6,16);
  l1->SetLineColor(kRed);
  // l1->Draw();

}// end of void coin_trig



