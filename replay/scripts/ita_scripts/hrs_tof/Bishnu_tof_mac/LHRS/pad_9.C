//////////////////////////////////////////////////////////////////
//07/11/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #9 
/////////////////////////////////////////////////////////////int rightwalk_second38()
//To copy a root file (coincidence experiment) to my directory

//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_9()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100592;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  
  /*
  // RF -s2 before calibration
  TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH1F *h9 = new TH1F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h9",pmt+1,pmt+49),cut_9);
  
  ////////////////////////////////////////////////
  //scintillator #9
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////
  
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.022,0.13);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,-0.02,0.12e-6,0.13);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.01,0.035);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.01,0.12e-6,0.035);
 l1p->SetLineColor(kRed);
 l1p->Draw();

 ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////

TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.022,0.13,300,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h9",pmt+1,pmt+49),cut_9,"colz");

 // For TProfile
  TProfile *h9t = new TProfile("h9t","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.022,0.13,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h9t",pmt+1,pmt+49),cut_9,"same");
 
  TF1 *ft = new TF1("ft","[0]+[1]*x",0.02,0.1);
  h9t->Fit("ft","R+");
  

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.04,0.17);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] ))>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,-0.04,0.12e-6,0.17);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.02,0.035);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] ))>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.02,0.12e-6,0.035);
 l1p->SetLineColor(kRed);
 l1p->Draw();
 
///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////
 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.02,0.035,300,0.1170e-6,0.122e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] )) :(L.tr.th[0])>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 // FOr TProfile
 TProfile *h9t = new TProfile("h9t","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,-0.02,0.035,0.1170e-6,0.122e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] )) :(L.tr.th[0])>>h9t",pmt+1,pmt+49),cut_9p,"same");
  TF1 *ft = new TF1("ft","[0]+[1]*x",0,0.02);
   h9t->Fit("ft","R+");
  
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
 
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.04,0.17);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  ))>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,-0.04,0.12e-6,0.17);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.017,0.035);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0] ))>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.017,0.12e-6,0.035);
 l1p->SetLineColor(kRed);
 l1p->Draw();


 
 
///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////

  
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form("(L.s2.la_c[9]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  ))>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form("(L.s2.ra_c[9]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0] ))>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,700);
 l1p->SetLineColor(kRed);
 l1p->Draw();
  

  ///////////////////////
  //Correcting for ADC left
  //
  //////////////////////////
  
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,700,300,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  )) :(L.s2.la_c[9])>>h9",pmt+1,pmt+49),cut_9,"colz");
 // FOr tprofile

 TProfile *h9t = new TProfile("h9t","LHRS pad#9 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,700,0.11750e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  )) :(L.s2.la_c[9])>>h9t",pmt+1,pmt+49),cut_9,"same");
 TF1 *ft = new TF1("ft","[0]+[1]*x",275,400);
   h9t->Fit("ft","R+");
  
////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form("(L.s2.la_c[9]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  -1.14584e-12*L.s2.la_c[9]  ))>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,700);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,700);

 T->Draw(Form("(L.s2.ra_c[9]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0] -1.14584e-12*L.s2.la_c[9] ))>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,700);
 l1p->SetLineColor(kRed);
 l1p->Draw();

  ////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////
  
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH2F *h9 = new TH2F("h9","LHRS pad#9 Y RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

 T->Draw(Form("(L.tr.y[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  -1.14584e-12*L.s2.la_c[9]  ))>>h9",pmt+1,pmt+49),cut_9,"colz");
 TLine *l1 = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1->SetLineColor(kRed);
  l1->Draw();

  // FOr X' 
TCut cut_9p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9p = new TCanvas("c9p","c9p", 600,600);
 
 TH2F *h9p = new TH2F("h9p","LHRS pad#9 Y'vs RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

 T->Draw(Form("(L.tr.ph[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0] -1.14584e-12*L.s2.la_c[9] ))>>h9p",pmt+1,pmt+49),cut_9p,"colz");
 TLine *l1p = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
 l1p->Draw();
*/
 ///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  //////////////////////////////
   /////////////////////
TCut cut_9 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==9    &&  L.s2.nthit==1");

 TCanvas *c9 = new TCanvas("c9","c9", 600,600);
 
 TH1F *h9 = new TH1F("h9","LHRS pad#9 Y RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  8.36477e-09*L.tr.x[0] -3.63706e-08*L.tr.th[0]  -1.14584e-12*L.s2.la_c[9]  ))>>h9",pmt+1,pmt+49),cut_9);
   


 return 0;
}
