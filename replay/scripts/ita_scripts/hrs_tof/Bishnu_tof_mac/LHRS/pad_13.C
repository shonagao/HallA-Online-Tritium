//////////////////////////////////////////////////////////////////
//07/11/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #13 
/////////////////////////////////////////////////////////////int rightwalk_second38()
//To copy a root file (coincidence experiment) to my directory

//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_13()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100350; i<100685;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  
  /*
  // RF -s2 before calibration
  TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH1F *h13 = new TH1F("h13","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h13",pmt+5,pmt+53),cut_13);
  
  ////////////////////////////////////////////////
  //scintillator #13
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////(L.s2.la_c[13])

TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.34,0.48);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h13",pmt+5,pmt+53),cut_13,"colz");
  
TLine *l1 = new TLine(0.12e-6,0.34,0.12e-6,0.48);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For x'

TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.05,0.099);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
TLine *l1p = new TLine(0.12e-6,0.05,0.12e-6,0.099);
 l1p->SetLineColor(kRed);
  l1p->Draw();


///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////
TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.34,0.48,300,0.1180e-6,0.123e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h13",pmt+5,pmt+53),cut_13,"colz");
  
 // FOr tProfile

 TProfile *h13t = new TProfile("h13t","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.34,0.48,0.1180e-6,0.123e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h13t",pmt+5,pmt+53),cut_13,"same");
  
TF1 *ft = new TF1("ft","[0]+[1]*x",0.38,0.44);
 h13t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.34,0.48);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0]   ))>>h13",pmt+5,pmt+53),cut_13,"colz");
  
TLine *l1 = new TLine(0.12e-6,0.34,0.12e-6,0.48);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For x'

TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.05,0.099);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0]  ))>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
TLine *l1p = new TLine(0.12e-6,0.05,0.12e-6,0.099);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////
TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.05,0.099,300,0.120e-6,0.124e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0]  )) :(L.tr.th[0])>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
 //For TProfile
 TProfile *h13t = new TProfile("h13t","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.05,0.099,0.120e-6,0.124e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0]  )) :(L.tr.th[0])>>h13t",pmt+5,pmt+53),cut_13p,"same");
  
TF1 *ft = new TF1("ft","[0]+[1]*x",0.064,0.081);
 h13t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
 TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.34,0.48);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]   ))>>h13",pmt+5,pmt+53),cut_13,"colz");
  
TLine *l1 = new TLine(0.12e-6,0.34,0.12e-6,0.48);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For x'

TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.05,0.099);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]  ))>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
TLine *l1p = new TLine(0.12e-6,0.05,0.12e-6,0.099);
 l1p->SetLineColor(kRed);
  l1p->Draw();


 // //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////

 TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#13 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,600);

 T->Draw(Form("(L.s2.la_c[13]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]   ))>>h13",pmt+5,pmt+53),cut_13,"colz");
  
TLine *l1 = new TLine(0.12e-6,150,0.12e-6,600);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For x'

TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#13 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,600);

 T->Draw(Form("(L.s2.ra_c[13]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]  ))>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
TLine *l1p = new TLine(0.12e-6,150,0.12e-6,600);
 l1p->SetLineColor(kRed);
 l1p->Draw();

//////////////////////
  //Correcting for ADC left
  //
  //////////////////////////
  
TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#13 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,150,600,300,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]   )) :(L.s2.la_c[13])>>h13",pmt+5,pmt+53),cut_13,"colz");
  
 //FOr Tprofile

 TProfile *h13t = new TProfile("h13t","LHRS pad#13 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,150,600,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]   )) :(L.s2.la_c[13])>>h13t",pmt+5,pmt+53),cut_13,"same");
  TF1 *ft = new TF1("ft","[0]+[1]*x",265,410);
 h13t->Fit("ft","R+");

////////////////////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 ////////////////////////////////
 
 
 TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH2F *h13 = new TH2F("h13","LHRS pad#13 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,600);

 T->Draw(Form("(L.s2.la_c[13]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]  -9.13463e-13*L.s2.la_c[13]  ))>>h13",pmt+5,pmt+53),cut_13,"colz");
 // -9.13462e-13*L.s2.la_c[13]
 //-1.06023e-12*L.s2.la_c[13]

TLine *l1 = new TLine(0.12e-6,150,0.12e-6,600);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For x'

TCut cut_13p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13p = new TCanvas("c13p","c13p", 600,600);
 
 TH2F *h13p = new TH2F("h13p","LHRS pad#13 ADC_R RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,150,600);

 T->Draw(Form("(L.s2.ra_c[13]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]  -9.13462e-13*L.s2.la_c[13]   ))>>h13p",pmt+5,pmt+53),cut_13p,"colz");
  
TLine *l1p = new TLine(0.12e-6,150,0.12e-6,600);
 l1p->SetLineColor(kRed);
 l1p->Draw();
//////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  //////////////////////////////////////////
  */
TCut cut_13 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==13    &&  L.s2.nthit==1");

 TCanvas *c13 = new TCanvas("c13","c13", 600,600);
 
 TH1F *h13 = new TH1F("h13","LHRS pad#13 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 6.22836e-09*L.tr.x[0] -3.07010e-08*L.tr.th[0]  -1.23475e-12*L.s2.la_c[13]  ))>>h13",pmt+5,pmt+53),cut_13);
 //-9.13462e-13*L.s2.la_c[13]
 //  -1.06023e-12*L.s2.la_c[13] 



 return 0;
}
