//////////////////////////////////////////////////////////////////
//07/10/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #2 
////////////////////////////////////////////////////////////
int pad_2()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100499;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));
 }
  Int_t pmt = 8;
  /*
  // RF -s2 with out any correction
TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH1F *h2 = new TH1F("h2","pad #2 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h2",pmt-6,pmt+42),cut_2, "colz");
  

////////////////////////////////////////////////
  //scintillator #2
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.61,-0.48);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h2",pmt-6,pmt+42),cut_2, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.61,0.12e-6,-0.48);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.122,-0.09);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h2p",pmt-6,pmt+42),cut_2p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.122,0.12e-6,-0.09);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 

 ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,-0.61,-0.48,200,0.1190e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h2",pmt-6,pmt+42),cut_2, "colz");

 // FOr TProfile

 TProfile *h2t = new TProfile("h2t","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,-0.61,-0.48,0.1190e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h2t",pmt-6,pmt+42),cut_2, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.578,-0.52);
  h2t->Fit("ft","R+");


/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////
  
TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.61,-0.48);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] ))>>h2",pmt-6,pmt+42),cut_2, "colz");
  TLine *l1 = new TLine(0.122e-6,-0.61,0.122e-6,-0.48);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.122,-0.09);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0]   ))>>h2p",pmt-6,pmt+42),cut_2p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.122,0.12e-6,-0.09);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 /////////////////////////////////////// 
  
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,-0.122,-0.09,200,0.1160e-6,0.1205e-6);
 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0]   )) :  (L.tr.th[0])>>h2p",pmt-6,pmt+42),cut_2p, "colz");

 // For Tprofile

 TProfile *h2t = new TProfile("h2t","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,-0.122,-0.09,0.1160e-6,0.1205e-6);
 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0]   )) :  (L.tr.th[0])>>h2t",pmt-6,pmt+42),cut_2p, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.114,-0.097);
  h2t->Fit("ft","R+");


///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.61,-0.48);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]  ))>>h2",pmt-6,pmt+42),cut_2, "colz");
  TLine *l1 = new TLine(0.122e-6,-0.61,0.122e-6,-0.48);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.122,-0.09);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]   ))>>h2p",pmt-6,pmt+42),cut_2p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.122,0.12e-6,-0.09);
 l1p->SetLineColor(kRed);
  l1p->Draw(); 
  
 ///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////// 

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 ADC _L RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.la_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]  ))>>h2",pmt-6,pmt+42),cut_2, "colz");
  TLine *l1 = new TLine(0.122e-6,180,0.122e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 ADC_R RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,800);

 T->Draw(Form(" (L.s2.ra_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]   ))>>h2p",pmt-6,pmt+42),cut_2p, "colz");
 TLine *l1p = new TLine(0.12e-6,170,0.12e-6,800);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  // FOR adc sum
TCut cut_2s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2s = new TCanvas("c2s","c2s", 600,600);
 
 TH2F *h2s = new TH2F("h2s","pad #2 ADC _sum RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,1100);

 T->Draw(Form(" (L.s2.la_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]  ))>>h2s",pmt-6,pmt+42),cut_2s, "colz");
  TLine *l1s = new TLine(0.122e-6,170,0.122e-6,1100);
 l1s->SetLineColor(kRed);
 l1s->Draw();

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 ADC _L RF -s2 time ;RF- s2 timein sec  ; ",200,180,800,200,0.12050e-6,0.126e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]  )) : (L.s2.la_c[2])>>h2",pmt-6,pmt+42),cut_2, "colz");


 //FOr TProfile

TProfile *h2t = new TProfile("h2t","pad #2 ADC _L RF -s2 time ;RF- s2 timein sec  ; ",200,180,800,0.12050e-6,0.126e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0]  )) : (L.s2.la_c[2])>>h2t",pmt-6,pmt+42),cut_2, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",240,475);
  h2t->Fit("ft","R+");


////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 ///////////



TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 ADC _L RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.la_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2]  ))>>h2",pmt-6,pmt+42),cut_2, "colz");
  TLine *l1 = new TLine(0.122e-6,180,0.122e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 ADC_R RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,800);

 T->Draw(Form(" (L.s2.ra_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2]    ))>>h2p",pmt-6,pmt+42),cut_2p, "colz");
 TLine *l1p = new TLine(0.12e-6,170,0.12e-6,800);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  // FOR adc sum
TCut cut_2s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2s = new TCanvas("c2s","c2s", 600,600);
 
 TH2F *h2s = new TH2F("h2s","pad #2 ADC _sum RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,1100);

 T->Draw(Form(" (L.s2.la_c[2]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2]   ))>>h2s",pmt-6,pmt+42),cut_2s, "colz");
  TLine *l1s = new TLine(0.122e-6,170,0.122e-6,1100);
 l1s->SetLineColor(kRed);
 l1s->Draw();


 ////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////

TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH2F *h2 = new TH2F("h2","pad #2 Y vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.05,0.042);

 T->Draw(Form(" (L.tr.y[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2]  ))>>h2",pmt-6,pmt+42),cut_2, "colz");
 TLine *l1 = new TLine(0.122e-6,-0.05,0.122e-6,0.042);
 l1->SetLineColor(kRed);
  l1->Draw(); 
  // FOr X' vs RF -s2 time
TCut cut_2p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2p = new TCanvas("c2p","c2p", 600,600);
 
 TH2F *h2p = new TH2F("h2p","pad #2 Y' vs  RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.04,0.05);

 T->Draw(Form(" (L.tr.ph[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 1.60686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2]    ))>>h2p",pmt-6,pmt+42),cut_2p, "colz");
 TLine *l1p = new TLine(0.12e-6,-0.04,0.12e-6,0.05);
 l1p->SetLineColor(kRed);
  l1p->Draw();


*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
TCut cut_2 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==2    &&  L.s2.nthit==1");

 TCanvas *c2 = new TCanvas("c2","c2", 600,600);
 
 TH1F *h2 = new TH1F("h2","pad #2 Y vs RF -s2 time ;RF- s2 timein sec  ; ",400,0.10e-6,0.15e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 2.39686e-09*L.tr.x[0] -3.32898e-08*L.tr.th[0] -1.24122e-12*L.s2.la_c[2] ))>>h2",pmt-6,pmt+42),cut_2);
 //  + 1.60686e-09*L.tr.x[0]
 return 0;
}
