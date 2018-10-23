//////////////////////////////////////////////////////////////////
//07/09/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #4 
/////////////////////////////////////////////////////////////int rightwalk_second38()
int pad_4()
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
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH1F *h4 = new TH1F("h4","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h4",pmt-4,pmt+44),cut_4, "colz");
  

////////////////////////////////////////////////
  //scintillator #4
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 X vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.48,-0.3);

 T->Draw(Form(" ( L.tr.x[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,-0.48,0.12e-6,-0.3);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.09,-0.045);

 T->Draw(Form(" ( L.tr.th[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,-0.09,0.12e-6,-0.045);
 l1p->SetLineColor(kRed);
  l1p->Draw();

 ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.48,-0.3,200,0.11850e-6,0.1228e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : ( L.tr.x[0]) >>h4",pmt-4,pmt+44),cut_4, "colz");
 // For TProfile

 TProfile *h4t = new TProfile("h4t","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.48,-0.3,0.11850e-6,0.1228e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : ( L.tr.x[0]) >>h4t",pmt-4,pmt+44),cut_4, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.44,-0.35);
  h4t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 X vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.48,-0.3);

 T->Draw(Form(" ( L.tr.x[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  ))>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,-0.48,0.12e-6,-0.3);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.09,-0.045);

 T->Draw(Form(" ( L.tr.th[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] ))>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,-0.09,0.12e-6,-0.045);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////
 
TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.09,-0.045,200,0.1180e-6,0.1225e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] )) : ( L.tr.th[0]) >>h4p",pmt-4,pmt+44),cut_4p, "colz");

 // FOr tProfile

 TProfile *h4t = new TProfile("h4t","pad #4 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.09,-0.045,0.1180e-6,0.1225e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] )) : ( L.tr.th[0]) >>h4t",pmt-4,pmt+44),cut_4p, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.08,-0.0555);
  h4t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////

 
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 X vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.48,-0.3);

 T->Draw(Form(" ( L.tr.x[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0]  ))>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,-0.48,0.12e-6,-0.3);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.09,-0.045);

 T->Draw(Form(" ( L.tr.th[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] -4.08767e-08*L.tr.th[0] ))>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,-0.09,0.12e-6,-0.045);
 l1p->SetLineColor(kRed);
  l1p->Draw();

  ///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////// 
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 ADC_L vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" ( L.s2.la_c[4]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0]  ))>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,1100);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For ADC right vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 ADC_R vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,800);

 T->Draw(Form(" ( L.s2.ra_c[4]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] -4.08767e-08*L.tr.th[0] ))>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,170,0.12e-6,800);
 l1p->SetLineColor(kRed);
 l1p->Draw();

 // FOr ADC sum

TCut cut_4s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4s = new TCanvas("c4s","c4s", 600,600);
 
 TH2F *h4s = new TH2F("h4s","pad #4 ADC_sum vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,450,1150);

 T->Draw(Form(" ( L.s2.la_c[4]) + L.s2.ra_c[4] :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0]  ))>>h4s",pmt-4,pmt+44),cut_4s, "colz");

 TLine *l1s = new TLine(0.12e-6,450,0.12e-6,1150);
 l1s->SetLineColor(kRed);
  l1s->Draw();

///////////////////////
  //Correcting for ADC left
  //
  ///////////////////////////
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 ADC_L vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,180,800,200,0.11950e-6,0.1235e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0]  )) :  ( L.s2.la_c[4]) >>h4",pmt-4,pmt+44),cut_4, "colz");

 // FOr TProfile
 TProfile *h4t = new TProfile("h4t","pad #4 ADC_L vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,180,800,0.11950e-6,0.1235e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0]  )) :  ( L.s2.la_c[4]) >>h4t",pmt-4,pmt+44),cut_4, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",250,525);
  h4t->Fit("ft","R+");


////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

/////////////////////

 
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 ADC_L vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" ( L.s2.la_c[4]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4]   ))>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,1100);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For ADC right vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 ADC_R vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,170,800);

 T->Draw(Form(" ( L.s2.ra_c[4]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4] ))>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,170,0.12e-6,800);
 l1p->SetLineColor(kRed);
 l1p->Draw();

 // FOr ADC sum

TCut cut_4s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4s = new TCanvas("c4s","c4s", 600,600);
 
 TH2F *h4s = new TH2F("h4s","pad #4 ADC_sum vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,450,1150);

 T->Draw(Form(" ( L.s2.la_c[4]) + L.s2.ra_c[4] :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4]  ))>>h4s",pmt-4,pmt+44),cut_4s, "colz");

 TLine *l1s = new TLine(0.12e-6,450,0.12e-6,1150);
 l1s->SetLineColor(kRed);
  l1s->Draw();

////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////

TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH2F *h4 = new TH2F("h4","pad #4 Y  vs  RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.05,0.042);

 T->Draw(Form(" ( L.tr.y[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4]   ))>>h4",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1 = new TLine(0.12e-6,-0.05,0.12e-6,0.042);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For ADC right vs RF -s2 time

TCut cut_4p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4p = new TCanvas("c4p","c4p", 600,600);
 
 TH2F *h4p = new TH2F("h4p","pad #4 Y'  vs  RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.04,0.05);

 T->Draw(Form(" ( L.tr.ph[0]) :((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.67284e-09*L.tr.x[0] -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4] ))>>h4p",pmt-4,pmt+44),cut_4, "colz");

 TLine *l1p = new TLine(0.12e-6,-0.04,0.12e-6,0.05);
 l1p->SetLineColor(kRed);
 l1p->Draw();
*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
TCut cut_4 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==4    &&  L.s2.nthit==1");

 TCanvas *c4 = new TCanvas("c4","c4", 600,600);
 
 TH1F *h4 = new TH1F("h4","pad #4 Y  vs  RF -s2 time ;RF- s2 timein sec  ; ",400,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.67284e-09*L.tr.x[0]  -4.08767e-08*L.tr.th[0] -1.16133e-12*L.s2.la_c[4]))>>h4",pmt-4,pmt+44),cut_4, "colz");



 return 0;
}
