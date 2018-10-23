//////////////////////////////////////////////////////////////////
//07/09/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #5 
/////////////////////////////////////////////////////////////int rightwalk_second38()
int pad_5()
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
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH1F *h5 = new TH1F("h5","pad #5 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h5",pmt-3,pmt+45),cut_5, "colz");
  
////////////////////////////////////////////////
  //scintillator #5
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////
 
  
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.4,-0.2);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h5",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.4,0.12e-6,-0.2);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_5p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5p = new TCanvas("c5p","c5p", 600,600);
 
 TH2F *h5p = new TH2F("h5p","pad #5 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.08,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h5p",pmt-3,pmt+45),cut_5p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.08,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw();


  ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////
  TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.38,-0.22,200,0.1186e-6,0.123e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h5",pmt-3,pmt+45),cut_5, "colz");

 // For Tprofile
TProfile *h5t = new TProfile("h5t","pad #5 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.38,-0.22,0.1186e-6,0.123e-6);

 T->Draw(Form("  ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) : (L.tr.x[0])>>h5t",pmt-3,pmt+45),cut_5, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.34,-0.26);
  h5t->Fit("ft","R+");


/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////



  
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.4,-0.2);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0]   ))>>h5",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.4,0.12e-6,-0.2);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_5p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5p = new TCanvas("c5p","c5p", 600,600);
 
 TH2F *h5p = new TH2F("h5p","pad #5 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.08,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0]  ))>>h5p",pmt-3,pmt+45),cut_5p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.08,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////
 
TCut cut_5p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5p = new TCanvas("c5p","c5p", 600,600);
 
 TH2F *h5p = new TH2F("h5p","pad #5 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.08,-0.03,200,0.1170e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0]  )) :  (L.tr.th[0])>>h5p",pmt-3,pmt+45),cut_5p, "colz");

 // For Tprofile

TProfile *h5t = new TProfile("h5t","pad #5 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.08,-0.03,0.1170e-6,0.122e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0]  )) :  (L.tr.th[0])>>h5t",pmt-3,pmt+45),cut_5p, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",-0.063,-0.043);
  h5t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////

 
  
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 X vs RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.4,-0.2);

 T->Draw(Form(" (L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    ))>>h5",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1 = new TLine(0.12e-6,-0.4,0.12e-6,-0.2);
 l1->SetLineColor(kRed);
  l1->Draw();
  // For X' vs RF -s2 time

TCut cut_5p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5p = new TCanvas("c5p","c5p", 600,600);
 
 TH2F *h5p = new TH2F("h5p","pad #5 X' vs  RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.08,-0.03);

 T->Draw(Form(" (L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0]  -4.07972e-08*L.tr.th[0]  ))>>h5p",pmt-3,pmt+45),cut_5p, "colz");
  TLine *l1p = new TLine(0.12e-6,-0.08,0.12e-6,-0.03);
 l1p->SetLineColor(kRed);
  l1p->Draw();

 ///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 ////////////////////////////////////////

TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.la_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    ))>>h5",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1 = new TLine(0.12e-6,180,0.12e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw();
  // FOr Right ADC
TCut cut_5r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5r = new TCanvas("c5r","c5r", 600,600);
 
 TH2F *h5r = new TH2F("h5r","pad #5 ADC_R vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.ra_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    ))>>h5r",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1r = new TLine(0.12e-6,180,0.12e-6,800);
 l1r->SetLineColor(kRed);
  l1r->Draw();

  // FOr ADC sum
TCut cut_5s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5s = new TCanvas("c5s","c5s", 600,600);
 
 TH2F *h5s = new TH2F("h5s","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,1100);

 T->Draw(Form(" (L.s2.la_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    ))>>h5s",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1s = new TLine(0.12e-6,1800,0.12e-6,1200);
 l1s->SetLineColor(kRed);
  l1s->Draw();

///////////////////////
  //Correcting for ADC left
  //
  ///////////////////////////

  
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,200,700,200,0.11850e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    )) : (L.s2.la_c[5])>>h5",pmt-3,pmt+45),cut_5, "colz");
 // FOr TProfile

 TProfile *h5t = new TProfile("h5t","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,200,700,0.11850e-6,0.1235e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]    )) : (L.s2.la_c[5])>>h5t",pmt-3,pmt+45),cut_5, "same");
TF1 *ft = new TF1("ft","[0]+[1]*x",260,500);
  h5t->Fit("ft","R+");

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.la_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]  -1.24632e-12*L.s2.la_c[5]))>>h5",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1 = new TLine(0.12e-6,180,0.12e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw();
  // FOr Right ADC
TCut cut_5r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5r = new TCanvas("c5r","c5r", 600,600);
 
 TH2F *h5r = new TH2F("h5r","pad #5 ADC_R vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form(" (L.s2.ra_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]  -1.24632e-12*L.s2.la_c[5]   ))>>h5r",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1r = new TLine(0.12e-6,180,0.12e-6,800);
 l1r->SetLineColor(kRed);
  l1r->Draw();

  // FOr ADC sum
TCut cut_5s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5s = new TCanvas("c5s","c5s", 600,600);
 
 TH2F *h5s = new TH2F("h5s","pad #5 ADC_L vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,1100);

 T->Draw(Form(" (L.s2.la_c[5]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0] -1.24632e-12*L.s2.la_c[5]))>>h5s",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1s = new TLine(0.12e-6,180,0.12e-6,1200);
 l1s->SetLineColor(kRed);
  l1s->Draw();

////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////

TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH2F *h5 = new TH2F("h5","pad #5 Y vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.05,0.042);

 T->Draw(Form(" (L.tr.y[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]  -1.24632e-12*L.s2.la_c[5]))>>h5",pmt-3,pmt+45),cut_5, "colz");
 TLine *l1 = new TLine(0.12e-6,-0.05,0.12e-6,0.042);
 l1->SetLineColor(kRed);
  l1->Draw();
  // FOr Right ADC
TCut cut_5r =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5r = new TCanvas("c5r","c5r", 600,600);
 
 TH2F *h5r = new TH2F("h5r","pad #5 Y' vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.04,0.05);

 T->Draw(Form(" (L.tr.ph[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]  -1.24632e-12*L.s2.la_c[5]   ))>>h5r",pmt-3,pmt+45),cut_5, "colz");
  TLine *l1r = new TLine(0.12e-6,-0.04,0.12e-6,0.05);
 l1r->SetLineColor(kRed);
  l1r->Draw();
*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
TCut cut_5 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==5    &&  L.s2.nthit==1");

 TCanvas *c5 = new TCanvas("c5","c5", 600,600);
 
 TH1F *h5 = new TH1F("h5","pad #5 Y vs RF -s2 time ;RF- s2 timein sec  ; ",400,0.10e-6,0.15e-6);

 T->Draw(Form("((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 +  7.88134e-09*L.tr.x[0] -4.07972e-08*L.tr.th[0]  -1.24632e-12*L.s2.la_c[5]))>>h5",pmt-3,pmt+45),cut_5);


  return 0;
}
