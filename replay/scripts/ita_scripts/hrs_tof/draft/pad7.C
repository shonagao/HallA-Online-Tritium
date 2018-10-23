// 09/26/2018
// Since RHRS timing is too much in comparision with LHRS
// trying to see RHRRS RF-S2 time for the second time.
// Following exactely same steps from LHRS ie pad_7.c 
// Author Bishnu
// Pad# 7


int pad7()
{

  
 
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100670; i<100685;i++) /// these run goes 100200 to 100685
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
 Int_t pmt = 8; //variable

 /* 
  //RF -S2 time before calibration
 TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH1F *h8 = new TH1F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6);

T->Draw(Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h8",pmt+15,pmt+47),cut_8);

 */

/*
// X vs RF -s2 uncorrected

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Draw(Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h8",pmt+15,pmt+47),cut_8,"colz");

  TLine *l1 = new TLine(0.222e-6,-0.25,0.222e-6,-0.02);
 l1->SetLineColor(kRed);
  l1->Draw();

*/  
  ///////////// corrections for X

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.22,-0.06,300,0.2206e-6,0.224e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(R.tr.x[0]) >>h8",pmt+15,pmt+47),cut_8,"colz");
 // Fit a Tprofile
 
 TProfile *ht8 = new TProfile("ht8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.22,-0.06,0.2206e-6,0.224e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(R.tr.x[0]) >>ht8",pmt+15,pmt+47),cut_8,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.2,-0.1);
  ht8->Fit("ft","R+");
  
  //// X vs RF -s2 (S2 is corrected for X)
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 nt corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Draw(Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.25,0.2235e-6,-0.02);
 l1->SetLineColor(kRed);
  l1->Draw();

  /*
 

  /// X' vs RF -s2 now s2 is already corrected for X

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.05,0);

 T->Draw(Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.05,0.2235e-6,0);
 l1->SetLineColor(kRed);
  l1->Draw();
  
  /// Correction for X'

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.05,0,300,0.2218e-6,0.2252e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  ) :(R.tr.th[0])>>h8",pmt+15,pmt+47),cut_8,"colz");
 // For TProfile
TProfile *ht8 = new TProfile("ht8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.05,0,0.2218e-6,0.2252e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  ) :(R.tr.th[0])>>ht8",pmt+15,pmt+47),cut_8,"same");


TF1 *ft = new TF1("ft","[0]+[1]*x",-0.037,-0.015);
  ht8->Fit("ft","R+");
  
  // Checking for X' vs RF -s2

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.05,0);

 T->Draw(Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  + 4.27527e-08*R.tr.th[0])>>h8",pmt+15,pmt+47),cut_8,"colz");
 //   p0           2.22662e-07 
  TLine *l1 = new TLine(0.2228e-6,-0.05,0.2228e-6,0);
 l1->SetLineColor(kRed);
  l1->Draw();

  
  // just check the X vs RF -s2

 TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Draw(Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.25,0.2235e-6,-0.02);
 l1->SetLineColor(kRed);
  l1->Draw();


  // ADCL vs RF -s2 

 TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,120,350);

 T->Draw(Form("(R.s2.la_c[7]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");
 
  TLine *l1 = new TLine(0.2235e-6,120,0.2235e-6,350);
 l1->SetLineColor(kRed);
  l1->Draw();
  // ADCL vs RF -S2

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,120,350,300,0.221e-6,0.2235e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] ):(R.s2.la_c[7]) >>h8",pmt+15,pmt+47),cut_8,"colz");
 
 // For tProfile
TProfile *ht8 = new TProfile("ht8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,120,350,0.221e-6,0.2235e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] ):(R.s2.la_c[7]) >>ht8",pmt+15,pmt+47),cut_8,"same");
 
TF1 *ft = new TF1("ft","[0]+[1]*x",185,250);
  ht8->Fit("ft","R+");
  
  // checking ADCL vs RF -s2 after ADC correction


TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,120,350);

 T->Draw(Form("(R.s2.ra_c[7]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] -2.37582e-13*R.s2.la_c[7]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.22499e-07 
  TLine *l1 = new TLine(0.2235e-6,120,0.2235e-6,350);
 l1->SetLineColor(kRed);
  l1->Draw();


  /// R.tr.y[0]  vs RF-s2
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.04,0.05);

 T->Draw(Form("(R.tr.ph[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] -2.37582e-13*R.s2.la_c[7]  )>>h8",pmt+15,pmt+47),cut_8,"colz");

 // This is final Histogram counts vs RF -s2 
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH1F *h8 = new TH1F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");

 //+2.37582e-13*R.s2.la_c[7]
  
*/

 return 0;
}
