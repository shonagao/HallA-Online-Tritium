//////////////////////////////////////////////////////////////////
//07/11/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #10 
/////////////////////////////////////////////////////////////int rightwalk_second38()
//To copy a root file (coincidence experiment) to my directory

//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_10()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100410;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  
  /*
  // RF -s2 before calibration
  TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH1F *h10 = new TH1F("h10","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h10",pmt+2,pmt+50),cut_10);
  
  ////////////////////////////////////////////////
  //scintillator #10
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////



TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.04,0.27);

 T->Draw(Form(" (L.tr.x[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h10",pmt+2,pmt+50),cut_10,"colz");
 TLine *l1 = new TLine(0.12e-6,0.04,0.12e-6,0.27);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.01,0.06);

 T->Draw(Form(" (L.tr.th[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.01,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.04,0.27,300,0.1170e-6,0.122e-6);
 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :  (L.tr.x[0])>>h10",pmt+2,pmt+50),cut_10,"colz");

 // For tprofile

TProfile *h10t = new TProfile("h10t","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.04,0.27,0.1170e-6,0.122e-6);
T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :  (L.tr.x[0])>>h10t",pmt+2,pmt+50),cut_10,"same");


TF1 *ft = new TF1("ft","[0]+[1]*x",0.11,0.2);
 h10t->Fit("ft","R+");

  /////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.04,0.27);

 T->Draw(Form(" (L.tr.x[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]   ))>>h10",pmt+2,pmt+50),cut_10,"colz");
 TLine *l1 = new TLine(0.12e-6,0.04,0.12e-6,0.27);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.01,0.06);

 T->Draw(Form(" (L.tr.th[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] ))>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.01,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
 l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////
 TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.01,0.06,300,0.1160e-6,0.1205e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] )) :(L.tr.th[0])>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 // For Tprofile
 
 TProfile *h10t = new TProfile("h10t","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.01,0.06,0.1160e-6,0.1205e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] )) :(L.tr.th[0])>>h10t",pmt+2,pmt+50),cut_10p,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",0.012,0.035);
 h10t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.04,0.27);

 T->Draw(Form(" (L.tr.x[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0]  ))>>h10",pmt+2,pmt+50),cut_10,"colz");
 TLine *l1 = new TLine(0.12e-6,0.04,0.12e-6,0.27);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.01,0.06);

 T->Draw(Form(" (L.tr.th[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] -3.39145e-08*L.tr.th[0] ))>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.01,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
 l1p->Draw();

///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form(" (L.s2.la_c[10]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0]  ))>>h10",pmt+2,pmt+50),cut_10,"colz");
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form(" (L.s2.ra_c[10]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] -3.39145e-08*L.tr.th[0] ))>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,700);
 l1p->SetLineColor(kRed);
 l1p->Draw();
 ///////////////////////
  //Correcting for ADC left
  //
  //////////////////////////
  
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,700,300,0.11650e-6,0.121e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0]  )) : (L.s2.la_c[10])>>h10",pmt+2,pmt+50),cut_10,"colz");

 // For Tprofile

 TProfile *h10t = new TProfile("h10t","LHRS pad#10 ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,700,0.11650e-6,0.121e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0]  )) : (L.s2.la_c[10])>>h10t",pmt+2,pmt+50),cut_10,"same");



TF1 *ft = new TF1("ft","[0]+[1]*x",275,415);
 h10t->Fit("ft","R+");

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 ///////////
  
TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form(" (L.s2.la_c[10]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0] -1.3259e-12*L.s2.la_c[10]   ))>>h10",pmt+2,pmt+50),cut_10,"colz");
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form(" (L.s2.ra_c[10]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +  8.74297e-09*L.tr.x[0] -3.39145e-08*L.tr.th[0] -1.3259e-12*L.s2.la_c[10] ))>>h10p",pmt+2,pmt+50),cut_10p,"colz");
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,700);
 l1p->SetLineColor(kRed);
 l1p->Draw();
 

 ////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////

TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH2F *h10 = new TH2F("h10","LHRS pad#10 Y  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);
T->Draw(Form(" (L.tr.y[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0] -1.3259e-12*L.s2.la_c[10]   ))>>h10",pmt+2,pmt+50),cut_10,"colz");
TLine *l1 = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_10p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10p = new TCanvas("c10p","c10p", 600,600);
 
 TH2F *h10p = new TH2F("h10p","LHRS pad#10 Y'vs  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);
 T->Draw(Form(" (L.tr.ph[0]):((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0] -1.3259e-12*L.s2.la_c[10]   ))>>h10",pmt+2,pmt+50),cut_10,"colz");
TLine *l1p = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
 l1p->Draw();
*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////

TCut cut_10 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==10    &&  L.s2.nthit==1");

 TCanvas *c10 = new TCanvas("c10","c10", 600,600);
 
 TH1F *h10 = new TH1F("h10","LHRS pad#10 Y  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);
 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.74297e-09*L.tr.x[0]  -3.39145e-08*L.tr.th[0] -1.3259e-12*L.s2.la_c[10]   ))>>h10",pmt+2,pmt+50),cut_10);

 return 0;
}
