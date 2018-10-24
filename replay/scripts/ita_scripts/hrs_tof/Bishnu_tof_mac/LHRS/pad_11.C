//////////////////////////////////////////////////////////////////
//07/11/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #11 
/////////////////////////////////////////////////////////////int rightwalk_second38()
//To copy a root file (coincidence experiment) to my directory

//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_11()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100028; i<100120;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  
  /*
  // RF -s2 before calibration
  TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH1F *h11 = new TH1F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h11",pmt+3,pmt+51),cut_11);
  
  ////////////////////////////////////////////////
  //scintillator #10
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////
  
 TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.15,0.35);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,0.15,0.12e-6,0.35);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 

  
TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.01,0.07);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.12e-6,0.01,0.12e-6,0.07);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////(L.tr.x[0])
  
 TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.15,0.35,300,0.1163e-6,0.1182e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0]) >>h11",pmt+3,pmt+51),cut_11,"colz");
  
 // For TProfile

TProfile *h11t = new TProfile("h11t","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.15,0.35,0.1163e-6,0.1182e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0]) >>h11t",pmt+3,pmt+51),cut_11,"same");
  
TF1 *ft = new TF1("ft","[0]+[1]*x",0.193,0.285);
 h11t->Fit("ft","R+");
 
 
 /////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////(L.tr.x[0])
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.15,0.35);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12  + 7.29128e-09*L.tr.x[0]  ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,0.15,0.12e-6,0.35);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 

  
TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS with updated + 7.29128e-09*L.tr.x ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.01,0.07);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  +7.29128e-09*L.tr.x[0]))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.12e-6,0.01,0.12e-6,0.07);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////

TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.01,0.07,300,0.1188e-6,0.1204e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 7.29128e-09*L.tr.x[0])) :(L.tr.th[0])  >>h11p",pmt+3,pmt+51),cut_11p,"colz");

 // FOr tProfile

 TProfile *h11t = new TProfile("h11t","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.01,0.07,0.11880e-6,0.1204e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 7.29128e-09*L.tr.x[0])) :(L.tr.th[0])  >>h11t",pmt+3,pmt+51),cut_11p,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",0.03,0.05);
 h11t->Fit("ft","R+");
  
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////

  TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
  
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.15,0.35);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12  + 7.29128e-09*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,0.15,0.12e-6,0.35);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 


TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.01,0.07);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 7.29128e-09*L.tr.x[0] -3.8785e-08*L.tr.th[0] ))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.1185e-6,0.01,0.1185e-6,0.07);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
  //////////////////////////
  //Making corrections for X for the second time
  ////////////////////////////////////
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.15,0.35,300,0.1134e-6,0.1158e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12  + 7.29128e-09*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0)) :(L.tr.x[0]) >>h11",pmt+3,pmt+51),cut_11,"colz");
  
 // For TProfile

 
 TProfile *h11t = new TProfile("h11t","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.15,0.35,0.1134e-6,0.1158e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12  + 7.29128e-09*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0)) :(L.tr.x[0]) >>h11t",pmt+3,pmt+51),cut_11,"same");
  
TF1 *ft = new TF1("ft","[0]+[1]*x",0.21,0.28);
 h11t->Fit("ft","R+");
 
/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.15,0.35);

 T->Draw(Form("(L.tr.x[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0]    ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,0.15,0.12e-6,0.35);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 
  

TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.01,0.07);

 T->Draw(Form("(L.tr.th[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.624608e-8*L.tr.x[0] -3.8785e-08*L.tr.th[0]  ))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.1225e-6,0.01,0.1225e-6,0.07);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
  
///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,650);

 T->Draw(Form("(L.s2.la_c[11]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,650);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 


TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,650);

 T->Draw(Form("(L.s2.ra_c[11]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.624608e-8*L.tr.x[0] -3.8785e-08*L.tr.th[0]  ))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,650);
 l1p->SetLineColor(kRed);
  l1p->Draw();

   
///////////////////////
  //Correcting for ADC left
  //
  //////////////////////////
  
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,650,300,0.1216e-6,0.1228e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0)) :(L.s2.la_c[11]) >>h11",pmt+3,pmt+51),cut_11,"colz");

   // FOr tProfile


 TProfile *h11t = new TProfile("h11t","LHRS pad#11ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,650,0.1216e-6,0.1228e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0]  ))/2.0)) :(L.s2.la_c[11]) >>h11t",pmt+3,pmt+51),cut_11,"same");

 TF1 *ft = new TF1("ft","[0]+[1]*x",260,375);
  h11t->Fit("ft","R+");

  
  
////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 ///////////
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,650);

 T->Draw(Form("(L.s2.la_c[11]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0] -1.63227e-12*L.s2.la_c[11] ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,180,0.12e-6,650);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 


TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,180,650);

 T->Draw(Form("(L.s2.ra_c[11]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.626608e-8*L.tr.x[0] -3.8785e-08*L.tr.th[0] -1.63227e-12*L.s2.la_c[11] ))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.12e-6,180,0.12e-6,650);
 l1p->SetLineColor(kRed);
 l1p->Draw();
  
  /////////////////////////////////
  ///Making corrections for DC right
  //////////////////////////////////////

TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,650,300,0.1305e-6,0.1325e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.546406e-8*L.tr.x[0] -3.05008e-08*L.tr.th[0] -2.4577e-12*L.s2.la_c[11] )) :(L.s2.ra_c[11]) >>h11p",pmt+3,pmt+51),cut_11p,"colz");

 // FOr TProfile

 TProfile *h11t = new TProfile("h11t","LHRS pad#11 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,180,650,0.1305e-6,0.1325e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -(((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.546406e-8*L.tr.x[0] -3.05008e-08*L.tr.th[0] -2.4577e-12*L.s2.la_c[11] )) :(L.s2.ra_c[11]) >>h11t",pmt+3,pmt+51),cut_11p,"same");
TF1 *ft = new TF1("ft","[0]+[1]*x",260,360);
 h11t->Fit("ft","R+");


 ////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////

TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH2F *h11 = new TH2F("h11","LHRS pad#11 Y vs   RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

 T->Draw(Form("(L.tr.y[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.624608e-8*L.tr.x[0]  -3.8785e-08*L.tr.th[0] -1.63227e-12*L.s2.la_c[11] ))/2.0))>>h11",pmt+3,pmt+51),cut_11,"colz");
  
 TLine *l1 = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1->SetLineColor(kRed);
  l1->Draw();



  // For X' 


TCut cut_11p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");

 TCanvas *c11p = new TCanvas("c11p","c11p", 600,600);
 
 TH2F *h11p = new TH2F("h11p","LHRS pad#11 Y' vs  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,-0.06,0.06);

 T->Draw(Form("(L.tr.ph[0]) : ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 -( ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 1.624608e-8*L.tr.x[0] -3.8785e-08*L.tr.th[0] -1.63227e-12*L.s2.la_c[11] ))>>h11p",pmt+3,pmt+51),cut_11p,"colz");
  
 TLine *l1p = new TLine(0.12e-6,-0.06,0.12e-6,0.06);
 l1p->SetLineColor(kRed);
 l1p->Draw();
  
  */
  
   ///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////
  
TCut cut_11 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==11    &&  L.s2.nthit==1");
 TCanvas *c11 = new TCanvas("c11","c11", 600,600);
 
 TH1F *h11 = new TH1F("h11","LHRS pad#11ADC _L  RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);

 // T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12   + 1.549e-8*L.tr.x[0]  -6.8785e-08*L.tr.th[0] -3.63227e-12*L.s2.la_c[11]  ))/2.0))>>h11",pmt+3,pmt+51),cut_11);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 7.745e-9*L.tr.x[0]  -3.43925e-08*L.tr.th[0] -1.816135e-12*L.s2.la_c[11]  ))>>h11",pmt+3,pmt+51),cut_11);
 

 //+1.86898e-13*L.s2.ra_c[11] 
 return 0;
}

