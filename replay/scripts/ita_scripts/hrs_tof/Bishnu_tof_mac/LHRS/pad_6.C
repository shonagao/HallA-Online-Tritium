//////////////////////////////////////////////////////////////////
//07/09/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #6 
/////////////////////////////////////////////////////////////int rightwalk_second38()
int pad_6()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100492;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));
 }
  Int_t pmt = 8; 
  /*
  /// RF -f2 before correction
TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH1F *h6 = new TH1F("h6","pad #6 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h6",pmt-2,pmt+46),cut_6);
  
 ////////////////////////////////////////////////
  //scintillator #6
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.320,-0.1);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h6",pmt-2,pmt+46),cut_6,"colz");
TLine *l1 = new TLine(0.12e-6,-0.32,0.12e-6,-0.1);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' vs RF -s2 time

TCut cut_6p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6p = new TCanvas("c6p","c6p", 600,600);
 
 TH2F *h6p = new TH2F("h6p","pad #6 X' vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.065,-0.015);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h6p",pmt-2,pmt+46),cut_6p,"colz");
TLine *l1p = new TLine(0.12e-6,-0.065,0.12e-6,-0.015);
 l1p->SetLineColor(kRed);
  l1p->Draw();

  ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////

TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.3,-0.13,200,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0]) >>h6",pmt-2,pmt+46),cut_6,"colz");

// For TProfile
 TProfile *h6t = new TProfile("h6t","pad #6 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.3,-0.13,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0]) >>h6t",pmt-2,pmt+46),cut_6,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.25,-0.17);
  h6t->Fit("ft","R+");
  

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.320,-0.1);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]   ))>>h6",pmt-2,pmt+46),cut_6,"colz");
TLine *l1 = new TLine(0.122e-6,-0.32,0.122e-6,-0.1);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' vs RF -s2 time

TCut cut_6p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6p = new TCanvas("c6p","c6p", 600,600);
 
 TH2F *h6p = new TH2F("h6p","pad #6 X' vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.065,-0.015);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  9.39930e-09*L.tr.x[0]  ))>>h6p",pmt-2,pmt+46),cut_6p,"colz");
TLine *l1p = new TLine(0.122e-6,-0.065,0.122e-6,-0.015);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////

 TCut cut_6p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6p = new TCanvas("c6p","c6p", 600,600);
 
 TH2F *h6p = new TH2F("h6p","pad #6 X' vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.06,-0.015,200,0.1275e-6,0.1325e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  9.39930e-09*L.tr.x[0]  )) : (L.tr.th[0]) >>h6p",pmt-2,pmt+46),cut_6p,"colz");

 // For TProfile
 TProfile *h6t = new TProfile("h6t","pad #6 X' vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.06,-0.015,0.1275e-6,0.1325e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  9.39930e-09*L.tr.x[0]  )) : (L.tr.th[0]) >>h6t",pmt-2,pmt+46),cut_6p,"same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.048,-0.028);
  h6t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.320,-0.1);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   ))>>h6",pmt-2,pmt+46),cut_6,"colz");
TLine *l1 = new TLine(0.1205e-6,-0.32,0.1205e-6,-0.1);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' vs RF -s2 time

TCut cut_6p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6p = new TCanvas("c6p","c6p", 600,600);
 
 TH2F *h6p = new TH2F("h6p","pad #6 X' vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.065,-0.015);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  9.39930e-09*L.tr.x[0] -3.63437e-08*L.tr.th[0]   ))>>h6p",pmt-2,pmt+46),cut_6p,"colz");
TLine *l1p = new TLine(0.1206e-6,-0.065,0.1206e-6,-0.015);
 l1p->SetLineColor(kRed);
  l1p->Draw();


///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////

TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 ADC_L vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,720);

 T->Draw(Form("(L.s2.la_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   ))>>h6",pmt-2,pmt+46),cut_6,"colz");
TLine *l1 = new TLine(0.1205e-6,180,0.1205e-6,720);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr ADC_right

TCut cut_6r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6r = new TCanvas("c6r","c6r", 600,600);
 
 TH2F *h6r = new TH2F("h6r","pad #6 ADC_R vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,720);

 T->Draw(Form("(L.s2.ra_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   ))>>h6r",pmt-2,pmt+46),cut_6r,"colz");
TLine *l1r = new TLine(0.1205e-6,180,0.1205e-6,720);
 l1r->SetLineColor(kRed);
 l1r->Draw();
 // For ADC sum
TCut cut_6s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6s = new TCanvas("c6s","c6s", 600,600);
 
 TH2F *h6s = new TH2F("h6s","pad #6 ADC_sum vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,400,1100);

 T->Draw(Form("(L.s2.la_c[6] + L.s2.ra_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   ))>>h6s",pmt-2,pmt+46),cut_6s,"colz");
TLine *l1s = new TLine(0.1205e-6,400,0.1205e-6,1100);
 l1s->SetLineColor(kRed);
  l1s->Draw();

///////////////////////
  //Correcting for ADC left
  //
  ///////////////////////////
TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 ADC_L vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,225,700,200,0.118e-6,0.123e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   )) : (L.s2.la_c[6])  >>h6",pmt-2,pmt+46),cut_6,"colz");

 // FOr TProfile
 TProfile *h6t = new TProfile("h6t","pad #6 ADC_L vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,225,700,0.118e-6,0.123e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]   )) : (L.s2.la_c[6])  >>h6t",pmt-2,pmt+46),cut_6,"same");
TF1 *ft = new TF1("ft","[0]+[1]*x",270,475);
  h6t->Fit("ft","R+");

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////



TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 ADC_L vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,720);

 T->Draw(Form("(L.s2.la_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6] ))>>h6",pmt-2,pmt+46),cut_6,"colz");
TLine *l1 = new TLine(0.1205e-6,180,0.1205e-6,720);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr ADC_right

TCut cut_6r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6r = new TCanvas("c6r","c6r", 600,600);
 
 TH2F *h6r = new TH2F("h6r","pad #6 ADC_R vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,720);

 T->Draw(Form("(L.s2.ra_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6]))>>h6r",pmt-2,pmt+46),cut_6r,"colz");
TLine *l1r = new TLine(0.1205e-6,180,0.1205e-6,720);
 l1r->SetLineColor(kRed);
 l1r->Draw();
 // For ADC sum
TCut cut_6s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6s = new TCanvas("c6s","c6s", 600,600);
 
 TH2F *h6s = new TH2F("h6s","pad #6 ADC_sum vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,400,1100);

 T->Draw(Form("(L.s2.la_c[6] + L.s2.ra_c[6]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6]))>>h6s",pmt-2,pmt+46),cut_6s,"colz");
TLine *l1s = new TLine(0.1205e-6,400,0.1205e-6,1100);
 l1s->SetLineColor(kRed);
  l1s->Draw();

////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////
TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH2F *h6 = new TH2F("h6","pad #6 Y  vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.05,0.042);

 T->Draw(Form("(L.tr.y[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6] ))>>h6",pmt-2,pmt+46),cut_6,"colz");
 TLine *l1 = new TLine(0.1205e-6,-0.05,0.1205e-6,0.042);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For Y'
TCut cut_6r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6r = new TCanvas("c6r","c6r", 600,600);
 
 TH2F *h6r = new TH2F("h6r","pad #6 Y'  vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.04,0.05);

 T->Draw(Form("(L.tr.ph[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6]))>>h6r",pmt-2,pmt+46),cut_6r,"colz");
TLine *l1r = new TLine(0.1205e-6,-0.04,0.1205e-6,0.05);
 l1r->SetLineColor(kRed);
 l1r->Draw();
*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
TCut cut_6 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==6    &&  L.s2.nthit==1");

 TCanvas *c6 = new TCanvas("c6","c6", 600,600);
 
 TH1F *h6 = new TH1F("h6","pad #6 counts  vs RF -s2 time(corrected s2) ;RF- s2 timein sec  ; ",500,0.10e-6,0.15e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  9.39930e-09*L.tr.x[0]  -3.63437e-08*L.tr.th[0]  -1.22617e-12*L.s2.la_c[6] ))>>h6",pmt-2,pmt+46),cut_6);





 return 0;
}
