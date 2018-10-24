//////////////////////////////////////////////////////////////////
//07/09/2018
//Author Bishnu Pandey
//This is the ldeft HRS ssystem s2 detector 
//calibrating scintillator of s2 detector(left HRS ) one by one and finally they all will be  combined together 
// scintillator #3 
/////////////////////////////////////////////////////////////int rightwalk_second38()
int pad_3()
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
TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH1F *h3 = new TH1F("h3","pad #3 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h3",pmt-5,pmt+43),cut_3, "colz");
  

////////////////////////////////////////////////
  //scintillator #3
  //This part  is the  X vs RF -s2 time 
  //note s2 is not corrected. WE are taking uncorrected s2 and correcting with respect to X,X' and ADC
  /////////////////////////////////////////////


TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.55,-0.38);

 T->Draw(Form("( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.55,0.12e-6,-0.38);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo rX' vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.11,-0.065);

 T->Draw(Form("( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0)>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.11,0.12e-6,-0.065);
 l1p->SetLineColor(kRed);
  l1p->Draw();

 ///////////////////////////////////////////
  // correcting for x.....RF -s2 time vs X
  //I will have (fx) after fitting a TProfile and this function will be added on  /// s2 time
  /////////////////////////////////////////////////////////
  

TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.55,-0.38,200,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :( L.tr.x[0])>>h3",pmt-5,pmt+43),cut_3, "colz");

 // FOr tProfile

 TProfile *h3t = new TProfile("h3t","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.55,-0.38,0.1180e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - ((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0) :( L.tr.x[0])>>h3t",pmt-5,pmt+43),cut_3, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.51,-0.42);
  h3t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X only
 //including the correction form X ie s2 = s2+ f(x)
 //////////////////////////////////////////////////////////////

 
TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.55,-0.38);

 T->Draw(Form("( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]    ))>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.55,0.12e-6,-0.38);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo rX' vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.11,-0.065);

 T->Draw(Form("( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] ))>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.11,0.12e-6,-0.065);
 l1p->SetLineColor(kRed);
  l1p->Draw();

///////////////////////////////////////
 // correction for X'
 //.RF -s2 time vs X'
 //I will have (fx') after fitting a TProfilr and this function will be added on s2 time again
 //
 ///////////////////////////////////////  
TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.11,-0.065,200,0.1170e-6,0.1215e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] )) :( L.tr.th[0])>>h3p",pmt-5,pmt+43),cut_3p, "colz");
 // For tprofile

TProfile *h3t = new TProfile("h3t","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,-0.11,-0.065,0.1170e-6,0.1215e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] )) :( L.tr.th[0])>>h3t",pmt-5,pmt+43),cut_3p, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",-0.096,-0.077);
  h3t->Fit("ft","R+");

/////////////////////////////////////////////
 //making two histograms X and X' vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X'ie s2 = s2+ f(x) +f(x')
 //////////////////////////////////////////////////////////////


TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.55,-0.38);

 T->Draw(Form("( L.tr.x[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]    ))>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.55,0.12e-6,-0.38);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo rX' vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4 RF -s2 time(uncorrected s2) ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.11,-0.065);

 T->Draw(Form("( L.tr.th[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] -4.14532e-08*L.tr.th[0] ))>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.11,0.12e-6,-0.065);
 l1p->SetLineColor(kRed);
  l1p->Draw();


 ///////////////////////////////////////////
 //making tthree histograms ADC left, ADC right and ADC sum vs RF -s2 time 
 //s is corrected for X  as well as X'
 //including the correction form X  and X' and X againie s2 = s2+ f(x) +f(x') 
 //////////////////////////////////////// 

TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 ADC left vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form("( L.s2.la_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]    ))>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,180,0.12e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo ADC_R vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4  ADC -right vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,150,800);

 T->Draw(Form("( L.s2.ra_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] -4.14532e-08*L.tr.th[0] ))>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,150,0.12e-6,800);
 l1p->SetLineColor(kRed);
 l1p->Draw();
 // FOr ADC sum
TCut cut_3s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3s = new TCanvas("c3s","c3s", 600,600);
 
 TH2F *h3s = new TH2F("h3s","pad #4 ADC sum vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,450,1150);

 T->Draw(Form("( L.s2.la_c[3]+L.s2.ra_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]    ))>>h3s",pmt-5,pmt+43),cut_3s, "colz");
   TLine *l1s = new TLine(0.12e-6,450,0.12e-6,1150);
 l1s->SetLineColor(kRed);
 l1s->Draw();
  
///////////////////////
  //Correcting for ADC left
  //
  ///////////////////////////
  
TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 ADC left vs RF -s2 time ;RF- s2 timein sec  ; ",200,180,800,200,0.1170e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]    )) :( L.s2.la_c[3])>>h3",pmt-5,pmt+43),cut_3, "colz");
 // FOr TProfile

 
 TProfile *h3t = new TProfile("h3t","pad #4 ADC left vs RF -s2 time ;RF- s2 timein sec  ; ",200,180,800,0.1170e-6,0.1225e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]    )) :( L.s2.la_c[3])>>h3t",pmt-5,pmt+43),cut_3, "same");

TF1 *ft = new TF1("ft","[0]+[1]*x",250,520);
  h3t->Fit("ft","R+");

////////////////////
 //s2 is now corrected for ADC L
 //Making three Histograms fro ADC L,R and sum vs RF -s2
 /////////////

 
TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 ADC left vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,180,800);

 T->Draw(Form("( L.s2.la_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]  -1.18518e-12*L.s2.la_c[3]   ))>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,180,0.12e-6,800);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo ADC_R vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4  ADC -right vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,150,800);

 T->Draw(Form("( L.s2.ra_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] -4.14532e-08*L.tr.th[0] -1.18518e-12*L.s2.la_c[3] ))>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,150,0.12e-6,800);
 l1p->SetLineColor(kRed);
 l1p->Draw();
 // FOr ADC sum
TCut cut_3s =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3s = new TCanvas("c3s","c3s", 600,600);
 
 TH2F *h3s = new TH2F("h3s","pad #4 ADC sum vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,450,1150);

 T->Draw(Form("( L.s2.la_c[3]+L.s2.ra_c[3]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0] -1.18518e-12*L.s2.la_c[3]   ))>>h3s",pmt-5,pmt+43),cut_3s, "colz");
   TLine *l1s = new TLine(0.12e-6,450,0.12e-6,1150);
 l1s->SetLineColor(kRed);
 l1s->Draw();

  ////////////////////////////
 //R.tr.y[0] vs RF -s2 and R.tr.ph[0] vs RF -s2
 //
 //
  /////////////////////
  TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH2F *h3 = new TH2F("h3","pad #4 Y  vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200, -0.05,0.042);

 T->Draw(Form("( L.tr.y[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]  -1.18518e-12*L.s2.la_c[3]   ))>>h3",pmt-5,pmt+43),cut_3, "colz");
   TLine *l1 = new TLine(0.12e-6,-0.05,0.12e-6,0.042);
 l1->SetLineColor(kRed);
  l1->Draw();
 // Fo ADC_R vs Rf- s2 time

TCut cut_3p =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3p = new TCanvas("c3p","c3p", 600,600);
 
 TH2F *h3p = new TH2F("h3p","pad #4  Y' vs RF -s2 time ;RF- s2 timein sec  ; ",200,0.10e-6,0.15e-6,200,-0.04,0.05);

 T->Draw(Form("( L.tr.ph[0]): ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0  + 5.50452e-09*L.tr.x[0] -4.14532e-08*L.tr.th[0] -1.18518e-12*L.s2.la_c[3] ))>>h3p",pmt-5,pmt+43),cut_3p, "colz");
   TLine *l1p = new TLine(0.12e-6,-0.04,0.12e-6,0.05);
 l1p->SetLineColor(kRed);
 l1p->Draw();

*/
///////////////////////////////////////////////
   //This is the final s2 = s2+ f(x)+f(x')+f(ADC_L)
  //
  //I am going to plot the no of counts vs RF - s2 time
  /////////////////////////////
 TCut cut_3 =("DR.evtypebits>>2&1 &&  L.s2.t_pads[0]==3    &&  L.s2.nthit==1");

 TCanvas *c3 = new TCanvas("c3","c3", 600,600);
 
 TH1F *h3 = new TH1F("h3","pad #3 counts  vs RF -s2 time ;RF- s2 timein sec  ; ",400,0.10e-6,0.15e-6);

 T->Draw(Form(" ((LTDC.F1FirstHit[37] -LTDC.F1FirstHit[47])*56.23e-12 - (((LTDC.F1FirstHit[30]-LTDC.F1FirstHit[%d])*56.23e-12 +((LTDC.F1FirstHit[37]-LTDC.F1FirstHit[%d])*56.23e-12))/2.0 + 5.50452e-09*L.tr.x[0]   -4.14532e-08*L.tr.th[0]  -1.18518e-12*L.s2.la_c[3]   ))>>h3",pmt-5,pmt+43),cut_3);

 return 0;
}
