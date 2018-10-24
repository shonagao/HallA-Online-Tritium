//////////////////////////////////////////////////////////////////
//07/10/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #1 
////////////////////////////////////////////////////////////
int pad_1()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100200; i<100651;i++)/// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));
 }
  Int_t pmt = 8;
  /*
  // RF -s2 with out any correction
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH1F *h1 = new TH1F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h1",pmt-7,pmt+41),cut_1, "colz");
  
 
////////////////////////////////////////////////
  //scintillator #1
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.68,-0.53);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.68,0.12e-6,-0.53);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  
  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.135,-0.11);

 T->Draw(Form("( L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.135,0.12e-6,-0.11);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 
  
 ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////
  
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.70,-0.52,200,0.12e-6,0.1219e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h1",pmt-7,pmt+41),cut_1, "colz");

 // FOr tProfile
 TProfile *h1t = new TProfile("h1t","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.70,-0.52,0.12e-6,0.1219e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h1t",pmt-7,pmt+41),cut_1, "same");

 TF1 *ft = new TF1("ft","[0]+[1]*x",-0.62,-0.59);
  h1t->Fit("ft","R+");
  
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.68,-0.57);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0] ))>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.68,0.12e-6,-0.57);
 l1->SetLineColor(kRed);
  l1->Draw(); 

  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.133,-0.113);

 T->Draw(Form("( L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0]  ))>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.133,0.12e-6,-0.113);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 
  
///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 /////////////////////////////////////// 
 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",500,-0.133,-0.113,1000,0.1180e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0]  )) :( L.tr.th[0]) >>h1p",pmt-7,pmt+41),cut_1p, "colz");
 // For TProfile
TProfile *h1t = new TProfile("h1t","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",500,-0.133,-0.113,0.1180e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0]  )) :( L.tr.th[0]) >>h1t",pmt-7,pmt+41),cut_1p, "same");
 TF1 *ft = new TF1("ft","[0]+[1]*x",-0.126,-0.118);
  h1t->Fit("ft","R+");

  
///////////////////////////////////////
 // making twp histograms X vs RF -s2 time and X' RF -s2 time
 //
 ///////////////////////////////////////

TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.68,-0.57);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] ))>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.68,0.12e-6,-0.57);
 l1->SetLineColor(kRed);
  l1->Draw(); 

  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.133,-0.113);

 T->Draw(Form("( L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0]  ))>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.133,0.12e-6,-0.113);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 


  //// ///////////////////////////////
  // making corrections for X for the second time
  ////////////////////////////////////
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,-0.68,-0.57,1000,0.1185e-6,0.1225e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] )) :(L.tr.x[0])>>h1",pmt-7,pmt+41),cut_1, "colz");

 // FOr TProfile
TProfile *h1t = new TProfile("h1t","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.68,-0.57,0.1185e-6,0.1225e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.04510e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] )) :(L.tr.x[0])>>h1t",pmt-7,pmt+41),cut_1, "same");

 TF1 *ft = new TF1("ft","[0]+[1]*x",-0.613,-0.59);
  h1t->Fit("ft","R+");
///////////////////////////////////////
 // making twp histograms X vs RF -s2 time and X' RF -s2 time
 // including corrections for X for the second time
 ///////////////////////////////////////


TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.68,-0.57);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] ))>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.68,0.12e-6,-0.57);
 l1->SetLineColor(kRed);
  l1->Draw(); 

  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,-0.133,-0.113);

 T->Draw(Form("( L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] ))>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.133,0.12e-6,-0.113);
 l1p->SetLineColor(kRed);
  l1p->Draw();





  
///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////// 
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 ADC_L  RF -s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,160,700);

 T->Draw(Form(" (L.s2.la_c[1]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] ))>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,160,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw(); 

  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1ADC_R RF -s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,160,700);

 T->Draw(Form("( L.s2.ra_c[1]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0]  ))>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,160,0.12e-6,700);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  

////////////////////
 //s2 is now correcting for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

 
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 ADC_L  RF -s2 time ;RF- s2 timein sec  ; ",500,200,500,200,0.12e-6,0.124e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] )) : (L.s2.la_c[1])>>h1",pmt-7,pmt+41),cut_1, "colz");

 // FOr TProfile
 
 TProfile *h1t = new TProfile("h1t","pad #1 ADC_L  RF -s2 time ;RF- s2 timein sec  ; ",400,200,500,0.12e-6,0.124e-6);
 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] )) : (L.s2.la_c[1])>>h1t",pmt-7,pmt+41),cut_1, "same");

 TF1 *ft = new TF1("ft","[0]+[1]*x",300,375);
  h1t->Fit("ft","R+");
  
  ///////////////////////////////
  //Two histograms after ADC correction
  ////////////////////////////////////

  
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH2F *h1 = new TH2F("h1","pad #1 ADC_L  RF -s2 time(-1.56e-12 updated)  ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,160,700);

 T->Draw(Form(" (L.s2.la_c[1]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0] -1.39426e-12 *L.s2.la_c[1]   ))>>h1",pmt-7,pmt+41),cut_1, "colz");
   TLine *l1 = new TLine(0.12e-6,160,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw(); 

  /// FOr X'vs Rf-s2 
TCut cut_1p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1p = new TCanvas("c1p","c1p", 600,600);
 
 TH2F *h1p = new TH2F("h1p","pad #1ADC_R RF -s2 time(-1.56e-12 ) ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,200,160,700);

 T->Draw(Form("( L.s2.ra_c[1]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -3.53783e-08*L.tr.th[0]  -1.39426e-12 *L.s2.la_c[1] ))>>h1p",pmt-7,pmt+41),cut_1p, "colz");
   TLine *l1p = new TLine(0.12e-6,160,0.12e-6,700);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  ///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
  
  */
TCut cut_1 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==1    &&  L.s2.nthit==1");

 TCanvas *c1 = new TCanvas("c1","c1", 600,600);
 
 TH1F *h1 = new TH1F("h1","pad #1 + 9.5982e-09*L.tr.x[0] -539783e-08*L.tr.th[0] - 1.592e-12 *L.s2.la_c[1]  )  ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);





   T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 9.5982e-09*L.tr.x[0] -5.3783e-08*L.tr.th[0] - 1.592e-12 *L.s2.la_c[1]))>>h1",pmt-7,pmt+41),cut_1);
  
   /// original 07/19/2018
   // 9.82e-09*L.tr.x[0] -4.9783e-08*L.tr.th[0] - 1.392e-12 *L.s2.la_c[1]))>>h1",pmt-7,pmt+41),cut_1);
 return 0;
}
// T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.44106e-09*L.tr.x[0] -5.8783e-08*L.tr.th[0] - 1.39426e-12 *L.s2.la_c[1]))>>h1",pmt-7,pmt+41),cut_1);
