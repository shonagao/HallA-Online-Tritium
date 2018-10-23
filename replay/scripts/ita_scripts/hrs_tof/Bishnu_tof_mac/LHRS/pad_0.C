//////////////////////////////////////////////////////////////////
//07/10/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #0 first scintillator
////////////////////////////////////////////////////////////
int pad_0()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100028; i<100386;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));
 }
  Int_t pmt = 8;
  /*
  // RF -s2 with out any correction
TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH1F *h0 = new TH1F("h0","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h0",pmt-8,pmt+40),cut_0, "colz");
  
  
////////////////////////////////////////////////
  //scintillator #0
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////
 
TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH2F *h0 = new TH2F("h0","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.75,-0.5);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h0",pmt-8,pmt+40),cut_0, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.75,0.12e-6,-0.5);
 l1->SetLineColor(kRed);
  l1->Draw();  
  // FOr X'
TCut cut_0p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0p = new TCanvas("c0p","c0p", 600,600);
 
 TH2F *h0p = new TH2F("h0p","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.19,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h0p",pmt-8,pmt+40),cut_0p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.19,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 
  
///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////
  TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH2F *h0 = new TH2F("h0","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",300,-0.75,-0.5,200,0.1170e-6,0.12e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0):(L.tr.x[0])>>h0",pmt-8,pmt+40),cut_0, "colz");

 //for Tprofile
TProfile *h0t = new TProfile("h0t","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",300,-0.75,-0.5,0.1170e-6,0.12e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0):(L.tr.x[0])>>h0t",pmt-8,pmt+40),cut_0, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.69,-0.66);
  h0t->Fit("ft","R+");
  
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH2F *h0 = new TH2F("h0","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.75,-0.5);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 2.07570e-8*L.tr.x[0]  ))>>h0",pmt-8,pmt+40),cut_0, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.75,0.12e-6,-0.5);
 l1->SetLineColor(kRed);
  l1->Draw();  
  // FOr X'
TCut cut_0p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0p = new TCanvas("c0p","c0p", 600,600);
 
 TH2F *h0p = new TH2F("h0p","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.19,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 2.07570e-08*L.tr.x[0] ))>>h0p",pmt-8,pmt+40),cut_0p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.19,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 
///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 /////////////////////////////////////// 

 
TCut cut_0p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0p = new TCanvas("c0p","c0p", 600,600);
 
 TH2F *h0p = new TH2F("h0p","pad #0 X' ;RF- s2 timein sec  ; ",300,-0.19,-0.03,200,0.1170e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 2.07570e-08*L.tr.x[0] )) : (L.tr.th[0])>>h0p",pmt-8,pmt+40),cut_0p, "colz");
 // FOr TProfile
 TProfile *h0t = new TProfile("h0t","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",300,-0.19,-0.03,0.1170e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 2.07570e-08*L.tr.x[0] )) : (L.tr.th[0])>>h0t",pmt-8,pmt+40),cut_0p, "same");
 
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.15,-0.14);
  h0t->Fit("ft","R+");
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s2 is corrected for X and X'
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////
  
TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH2F *h0 = new TH2F("h0","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.75,-0.5);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.92399e-09*L.tr.x[0] -3.51622e-08*L.tr.th[0] ))>>h0",pmt-8,pmt+40),cut_0, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.75,0.12e-6,-0.5);
 l1->SetLineColor(kRed);
  l1->Draw();  
  // FOr X'
TCut cut_0p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0p = new TCanvas("c0p","c0p", 600,600);
 
 TH2F *h0p = new TH2F("h0p","pad #0 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,-0.19,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.92399e-09*L.tr.x[0] -3.51622e-08*L.tr.th[0]  ))>>h0p",pmt-8,pmt+40),cut_0p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.19,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw();

  

////////////////////
 //s2 is now correcting for  ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH2F *h0 = new TH2F("h0","pad #0 ADC_L RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,300,200,500);

 T->Draw(Form(" (L.s2.la_c[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 2.07570e-08*L.tr.x[0] -3.51622e-08*L.tr.th[0] ))>>h0",pmt-8,pmt+40),cut_0, "colz");
 

 // -1.69994e-08 

 */
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////

TCut cut_0 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==0    &&  L.s2.nthit==1");

 TCanvas *c0 = new TCanvas("c0","c0", 600,600);
 
 TH1F *h0 = new TH1F("h0","pad#0  with L.t.x and L.tr.th   ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 // T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.92399*L.tr.x[0] -4.1622e-08*L.tr.th[0] -2.3563e-12*L.s2.la_c[0] ))>>h0",pmt-8,pmt+40),cut_0, "colz");
 // + 1.92399e-09*L.tr.x[0] -3.51622e-08*L.tr.th[0] -1.3563e-12*L.s2.la_c[0] 

  T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 2.0757e-8*L.tr.x[0] -2.69994e-08*L.tr.th[0] -1.4063e-12*L.s2.la_c[0]  ))>>h0",pmt-8,pmt+40),cut_0, "colz");



return 0;
}
