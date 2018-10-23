//////////////////////////////////////////////////////////////////
//07/05/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #8 
/////////////////////////////////////////////////////////////int rightwalk_second38()


//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_8()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100422;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  /*

  // RF -s2 before calibration
  TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH1F *h8 = new TH1F("h8","LHRS pad#8 RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h8",pmt,pmt+48),cut_8);
  
  ////////////////////////////////////////////////
  //scintillator #8
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," X vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.15,0.1);

  T->Draw(Form(" ( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h8",pmt,pmt+48),cut_8,"colz");

  TLine *l1 = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1->SetLineColor(kRed);
  l1->Draw();
  //// for   X' vs RF -s2 

 TCut cut_8p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8p = new TCanvas("c8p","c8p", 600,600);
 TH2F *h8p = new TH2F("h8p","X'vs  RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.035,0.025);

  T->Draw(Form(" ( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h8p",pmt,pmt+48),cut_8p,"colz");

  TLine *l1p = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
   ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////

 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8","  RF -s2 time vs X ;RF- s2 timein sec  ; ",200,-0.12,0.06,300,0.118e-6,0.122e-6);

  T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : ( L.tr.x[0]) >>h8",pmt,pmt+48),cut_8,"colz");
  // For TProfile
 TProfile *h8t = new TProfile("h8t","  RF -s2 time vs X ;RF- s2 timein sec  ; ",200,-0.12,0.06,0.118e-6,0.122e-6);

  T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : ( L.tr.x[0]) >>h8t",pmt,pmt+48),cut_8,"same");

  TF1 *ft = new TF1("ft","[0]+[1]*x",-0.07,0.01);
  h8t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," X vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.15,0.1);

  T->Draw(Form(" ( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0]  ))>>h8",pmt,pmt+48),cut_8,"colz");

  TLine *l1 = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1->SetLineColor(kRed);
  l1->Draw();
  //// for   X' vs RF -s2 

 TCut cut_8p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8p = new TCanvas("c8p","c8p", 600,600);
 TH2F *h8p = new TH2F("h8p","X'vs  RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.035,0.025);

  T->Draw(Form(" ( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.13305e-09*L.tr.x[0] ))>>h8p",pmt,pmt+48),cut_8p,"colz");

  TLine *l1p = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////

TCut cut_8p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8p = new TCanvas("c8p","c8p", 600,600);
 TH2F *h8p = new TH2F("h8p","X'vs  RF -s2 time ;RF- s2 timein sec  ; ",200,-0.025,0.014,300,0.126e-6,0.1308e-6);

  T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.13305e-09*L.tr.x[0] )) : ( L.tr.th[0])>>h8p",pmt,pmt+48),cut_8p,"colz");

  // FOr TProfile
TProfile *h8t = new TProfile("h8t","X'vs  RF -s2 time ;RF- s2 timein sec  ; ",200,-0.025,0.014,0.126e-6,0.1308e-6);

  T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.13305e-09*L.tr.x[0] )) : ( L.tr.th[0])>>h8t",pmt,pmt+48),cut_8p,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.015,0.002);
  h8t->Fit("ft","R+");


/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////

 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," X vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.15,0.1);

  T->Draw(Form(" ( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  ))>>h8",pmt,pmt+48),cut_8,"colz");

  TLine *l1 = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1->SetLineColor(kRed);
  l1->Draw();
  //// for   X' vs RF -s2 

 TCut cut_8p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8p = new TCanvas("c8p","c8p", 600,600);
 TH2F *h8p = new TH2F("h8p","X'vs  RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.035,0.025);

  T->Draw(Form(" ( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] ))>>h8p",pmt,pmt+48),cut_8p,"colz");

  TLine *l1p = new TLine(0.12e-6,-0.15,0.12e-6,0.1);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  ///////////////////////////////////
  // no need to do correction for X or X' they looks fine.
  ///////////////////////////////////////////

  
///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////

 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,650);

  T->Draw(Form(" ( L.s2.la_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  ))>>h8",pmt,pmt+48),cut_8,"colz");
TLine *l1 = new TLine(0.12e-6,150,0.12e-6,650);
 l1->SetLineColor(kRed);
  l1->Draw();
  /// FOr ADC right

TCut cut_8r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8r = new TCanvas("c8r","c8r", 600,600);
 TH2F *h8r = new TH2F("h8r"," ADC_R vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,650);

  T->Draw(Form(" ( L.s2.ra_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  ))>>h8r",pmt,pmt+48),cut_8r,"colz");
TLine *l1r = new TLine(0.12e-6,150,0.12e-6,650);
 l1r->SetLineColor(kRed);
  l1r->Draw();
  // For ADC sum

 TCut cut_8s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8s = new TCanvas("c8s","c8s", 600,600);
 TH2F *h8s = new TH2F("h8s"," ADC_sum vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,400,1100);

  T->Draw(Form(" ( L.s2.la_c[8]+L.s2.ra_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  ))>>h8s",pmt,pmt+48),cut_8,"colz");
TLine *l1s = new TLine(0.12e-6,400,0.12e-6,1100);
 l1s->SetLineColor(kRed);
  l1s->Draw();

 ///////////////////////
  //Correcting for ADC left
  //
  ///////////////////////////
 TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,175,550,300,0.13e-6,0.135e-6);

  T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  )) :( L.s2.la_c[8])>>h8",pmt,pmt+48),cut_8,"colz");
  // FOr TProfile
 TProfile *h8t = new TProfile("h8t"," ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,175,550,0.13e-6,0.135e-6);

  T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0]  )) :( L.s2.la_c[8])>>h8t",pmt,pmt+48),cut_8,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",275,430);
 h8t->Fit("ft","R+");

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,650);

  T->Draw(Form(" ( L.s2.la_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8",pmt,pmt+48),cut_8,"colz");
TLine *l1 = new TLine(0.1202e-6,150,0.1202e-6,650);
 l1->SetLineColor(kRed);
  l1->Draw();
  /// FOr ADC right

TCut cut_8r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8r = new TCanvas("c8r","c8r", 600,600);
 TH2F *h8r = new TH2F("h8r"," ADC_R vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,650);

  T->Draw(Form(" ( L.s2.ra_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8r",pmt,pmt+48),cut_8r,"colz");
TLine *l1r = new TLine(0.12e-6,150,0.12e-6,650);
 l1r->SetLineColor(kRed);
  l1r->Draw();
  // For ADC sum

 TCut cut_8s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8s = new TCanvas("c8s","c8s", 600,600);
 TH2F *h8s = new TH2F("h8s"," ADC_sum vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,400,1100);

  T->Draw(Form(" ( L.s2.la_c[8]+L.s2.ra_c[8]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8s",pmt,pmt+48),cut_8,"colz");
TLine *l1s = new TLine(0.12e-6,400,0.12e-6,1100);
 l1s->SetLineColor(kRed);
  l1s->Draw();
  
////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////
  
TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH2F *h8 = new TH2F("h8"," Y  vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

  T->Draw(Form(" ( L.tr.y[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8",pmt,pmt+48),cut_8,"colz");
TLine *l1 = new TLine(0.1204e-6,-0.06,0.1204e-6,0.06);
 l1->SetLineColor(kRed);
  l1->Draw();
  /// FOr ADC right

TCut cut_8r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8r = new TCanvas("c8r","c8r", 600,600);
 TH2F *h8r = new TH2F("h8r"," Y' vs RF -s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

  T->Draw(Form(" ( L.tr.ph[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8r",pmt,pmt+48),cut_8r,"colz");
TLine *l1r = new TLine(0.1204e-6,-0.06,0.1204e-6,0.06);
 l1r->SetLineColor(kRed);
  l1r->Draw();
*/
   ///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  //////////////////////////////
   /////////////////////

  
TCut cut_8 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==8    &&  L.s2.nthit==1");
 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 TH1F *h8 = new TH1F("h8"," counts  RF -s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);

  T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.13305e-09*L.tr.x[0] -3.54952e-08*L.tr.th[0] -6.57058e-13*L.s2.la_c[8]))>>h8",pmt,pmt+48),cut_8);


  

  return 0;
}
