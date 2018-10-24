//////////////////////////////////////////////////////////////////
//07/11/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all are combined 
// scintillator #12 
/////////////////////////////////////////////////////////////int rightwalk_second38()
//To copy a root file (coincidence experiment) to my directory

//[a-onl@aonl1 Rootfiles]$  cp /volatile/halla/triton/eep_Rootfiles/pass1/tritium_100429*.root ./

int pad_12()
{

  
  // gStyle->SetOptFit(1111110);
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100400; i<100470;i++) /// these run goes 90854 to 90862
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
  Int_t pmt = 8; 

  
  /*
  // RF -s2 before calibration
  TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH1F *h12 = new TH1F("h12","LHRS pad#11 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6);

T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h12",pmt+4,pmt+52),cut_12);
  
  ////////////////////////////////////////////////
  //scintillator #10
  //This part (2a and 2b) is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////

 TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.25,0.40);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h12",pmt+4,pmt+52),cut_12, "colz");
  TLine *l1 = new TLine(0.12e-6,0.25,0.12e-6,0.40);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.10e-6,0.15e-6,200,0.03,0.08);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h12p",pmt+4,pmt+52),cut_12p,"colz");
  TLine *l1p = new TLine(0.12e-6,0.03,0.12e-6,0.08);
 l1p->SetLineColor(kRed);
  l1p->Draw();
///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  //////////////////////////////////////////////////////////
  
TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.25,0.40,500,0.1190e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h12",pmt+4,pmt+52),cut_12, "colz");

 //For Tprofile
 TProfile *h12t = new TProfile("h12t","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",200,0.25,0.40,0.1190e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :(L.tr.x[0])>>h12t",pmt+4,pmt+52),cut_12, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",0.29,0.36);
 h12t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////


  
 TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,0.25,0.40);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.64972e-09*L.tr.x[0] ))>>h12",pmt+4,pmt+52),cut_12, "colz");
  TLine *l1 = new TLine(0.12e-6,0.25,0.12e-6,0.40);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,0.03,0.08);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.64972e-09*L.tr.x[0] ))>>h12p",pmt+4,pmt+52),cut_12p,"colz");
  TLine *l1p = new TLine(0.12e-6,0.03,0.12e-6,0.08);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.03,0.08,1000,0.1160e-6,0.121e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.64972e-09*L.tr.x[0] )) : (L.tr.th[0])>>h12p",pmt+4,pmt+52),cut_12p,"colz");

 // For Tprofile
TProfile *h12t = new TProfile("h12t","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,0.03,0.08,0.1160e-6,0.121e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 8.64972e-09*L.tr.x[0] )) : (L.tr.th[0])>>h12t",pmt+4,pmt+52),cut_12p,"same");

TF1 *ft = new TF1("ft","[0]+[1]*x",0.045,0.062);
 h12t->Fit("ft","R+");


/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////
 

TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,0.25,0.40);

 T->Draw(Form("(L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  ))>>h12",pmt+4,pmt+52),cut_12, "colz");
  TLine *l1 = new TLine(0.12e-6,0.25,0.12e-6,0.40);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,0.03,0.08);

 T->Draw(Form("(L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  ))>>h12p",pmt+4,pmt+52),cut_12p,"colz");
  TLine *l1p = new TLine(0.12e-6,0.03,0.12e-6,0.08);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
  // //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////////////////////////////
TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,150,600);

 T->Draw(Form("(L.s2.la_c[12]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  ))>>h12",pmt+4,pmt+52),cut_12, "colz");
  TLine *l1 = new TLine(0.12e-6,150,0.12e-6,600);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,150,600);

 T->Draw(Form("(L.s2.ra_c[12]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  ))>>h12p",pmt+4,pmt+52),cut_12p,"colz");
  TLine *l1p = new TLine(0.12e-6,150,0.12e-6,600);
 l1p->SetLineColor(kRed);
  l1p->Draw();
  
///////////////////////
  //Correcting for ADC left
  //
  //////////////////////////
TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,150,600,1000,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  )) :(L.s2.la_c[12])>>h12",pmt+4,pmt+52),cut_12, "colz");
 // For TProfile
TProfile *h12t = new TProfile("h12t","LHRS pad#12 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",300,150,600,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  )) :(L.s2.la_c[12])>>h12t",pmt+4,pmt+52),cut_12, "same");
 TF1 *ft = new TF1("ft","[0]+[1]*x",275,400);
 h12t->Fit("ft","R+");

////////////////////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 ////////////////////////////////

TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH2F *h12 = new TH2F("h12","LHRS pad#12 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,150,600);

 T->Draw(Form("(L.s2.la_c[12]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0] -1.14396e-12*L.s2.la_c[12]   ))>>h12",pmt+4,pmt+52),cut_12, "colz");
  TLine *l1 = new TLine(0.12e-6,150,0.12e-6,600);
 l1->SetLineColor(kRed);
  l1->Draw();

  // For X'

TCut cut_12p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12p = new TCanvas("c12p","c12p", 600,600);
 
 TH2F *h12p = new TH2F("h12p","LHRS pad#12 RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6,500,150,600);

 T->Draw(Form("(L.s2.ra_c[12]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  -1.14396e-12*L.s2.la_c[12]  ))>>h12p",pmt+4,pmt+52),cut_12p,"colz");
  TLine *l1p = new TLine(0.12e-6,150,0.12e-6,600);
 l1p->SetLineColor(kRed);
  l1p->Draw();





 */

///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////


TCut cut_12 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==12    &&  L.s2.nthit==1");

 TCanvas *c12 = new TCanvas("c12","c12", 600,600);
 
 TH1F *h12 = new TH1F("h12","LHRS pad#12 ADC_L RF -s2 time(uncorrected RF-s2 time ;RF- s2 timein sec  ; ",1000,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.94972e-09*L.tr.x[0] -3.05521e-08*L.tr.th[0]  -1.15547e-12*L.s2.la_c[12]  ))>>h12",pmt+4,pmt+52),cut_12, "colz");
 //+ 8.64972e-09*L.tr.x[0]


 return 0;
}
