// 09/26/2018
// Since RHRS timing is too much in comparision with LHRS
// trying to see RHRRS RF-S2 time for the second time.
// Following exactely same steps from LHRS ie pad_7.c 
// Author Bishnu 
// Modify Itabashi 09/27/2018
// Pad# 7

//TCut cut_tof(int i);

void  hrs_tof_correct ()
{

  
 
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100670; i<100680;i++) /// these run goes 100200 to 100685
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
 Int_t pmtl = 16; //PMT Left variable
 Int_t pmtr = 48;  // PMT Right variabl
 int ch; //channel
 ch=7; //select channel
 int chmax=16;
 //====== Defined Function ================//
 double p0[chmax];
 double p1[chmax];
 TCanvas *c[chmax];
 TCut cut[chmax];
 TH1F *htof[chmax];
 //=== X correction ====//
 TH2F *hx[chmax];
 TProfile *tpx[chmax];
 TF1 *fx[chmax];
 TH2F *hx_c[chmax];
 TH1F *htof_xc[chmax];
 //== X' correction =====//
 TCanvas *cxp[chmax];
 TH2F *hxp[chmax];
 TProfile *tpxp[chmax];
 TF1 *fxp[chmax];
 TH2F *hxp_c[chmax];
 TH1F *htof_xpc[chmax];
 //=== Y correction ====//
 TH2F *hy[chmax];
 TProfile *tpy[chmax];
 TF1 *fy[chmax];
 TH2F *hy_c[chmax];
 TH1F *htof_yc[chmax];
 //== Y' correction =====//
 TCanvas *cyp[chmax];
 TH2F *hyp[chmax];
 TProfile *tpyp[chmax];
 TF1 *fyp[chmax];
 TH2F *hyp_c[chmax];
 TH1F *htof_ypc[chmax];

 //===================== Cut condition===========================//

TCut cut_tof[chmax];
double tof_min[chmax],tof_max[chmax];

 tof_min[7]=0.225;
 tof_max[7]=0.227;


 cut_tof[7]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>(%lf)*1.0e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<(%lf)*1.0e-6 ",pmtl+7,pmtr+7,tof_min[7],pmtl+7,pmtr+7,tof_max[7]);



 //==============================================================//


  //RF -S2 time before calibration


 for(int i=0;i<chmax;i++){

   if(i==ch){


   c[i] = new TCanvas(Form("c%d",i),Form("c%d",i));
   c[i]->Divide(2,2);


 cut[i]=(Form("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==%d    &&  R.s2.nthit==1",i));
 


 //-------- TOF hist w/o calibration -----------//
 htof[i] = new TH1F(Form("htof[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6);
T->Project(Form("htof[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",pmtl+i,pmtr+i),cut[i]);

c[i]->cd(1);
 htof[i]->Draw();
  
//-------- TOF vs X hist --------------------------//
 hx[i] = new TH2F(Form("hx[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.25,-0.02);

T->Project(Form("hx[%d]",i),Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",pmtl+i,pmtr+i),cut[i]);
c[i]->cd(2);
 hx[i]->Draw();
//---------TProfile  Correction for X ----------------------//
 c[i]->cd(3);
 tpx[i] = new TProfile(Form("tpx[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 xposi ; ",i),200,-0.22,-0.06,tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);

 T->Project(Form("tpx[%d]",i),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(R.tr.x[0])",pmtl+i,pmtr+i),cut[i] &&cut_tof[i]);
 tpx[i]->Draw();

 //------ Fit TProfile Correction for X -------------------//
 //fx[i]= new TF1(Form("fx[%d]",i),"[0]+[1]*x",0.2206e-6,0.224e-6);
  fx[i]= new TF1(Form("fx[%d]",i),"[0]+[1]*x",tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 tpx[i]->Fit(Form("fx[%d]",i),"","");
 tpx[i]->GetYaxis()->SetRangeUser(tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 fx[i]->Draw("same"); 
  c[i]->cd(4);

 //----- Corrected for X -------------//  
 
 
 p0[i]=fx[i]->GetParameter(0);
 p1[i]=fx[i]->GetParameter(1);
 cout<<"Get parameter [0] is "<<p0[i]<<endl;
 cout<<"Get parameter [1] is "<<p1[i]<<endl;
 p1[i]=p1[i]*1.0e9; // [ns]


 hx_c[i] = new TH2F(Form("hx_c[%d]",i),Form("pad #%d RF -s2 time(s2 after corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Project(Form("hx_c[%d]",i),Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -R.tr.x[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i]);
  hx_c[i]->Draw();



 //--------- TOF affter corrected X ----------//

 htof_xc[i] = new TH1F(Form("htof_xc[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6);

 T->Project(Form("htof_xc[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -R.tr.x[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i] && cut_tof[i]);

 c[i]->cd(1);
 htof_xc[i]->SetLineColor(2);
 htof_xc[i]->Draw("same");

 //========== X' Correction ====================//

cxp[i] = new TCanvas(Form("cxp[%d]",i),Form("cxp[%d]",i));
   cxp[i]->Divide(2,2);


//-------- TOF vs X' hist --------------------------//
 hxp[i] = new TH2F(Form("hxp[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.05,0.0);

T->Project(Form("hxp[%d]",i),Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",pmtl+i,pmtr+i),cut[i]);
cxp[i]->cd(1);
 hxp[i]->Draw();


//---------TProfile  Correction for X' ----------------------//
 cxp[i]->cd(2);

 tpxp[i] = new TProfile(Form("tpxp[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;s2 x' position  ; ",i),200,-0.05,0.00,0.2*1.0e-6,0.25*1.0e-6);

 T->Project(Form("tpxp[%d]",i),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(R.tr.th[0])",pmtl+i,pmtr+i),cut[i]);
 tpxp[i]->Draw();
 //------ Fit TProfile Correction for X' -------------------//

 fxp[i]= new TF1(Form("fxp[%d]",i),"[0]+[1]*x",0.2206e-6,0.224e-6);
 tpxp[i]->Fit(Form("fxp[%d]",i),"","");
 p0[i]=fxp[i]->GetParameter(0);
 p1[i]=fxp[i]->GetParameter(1);

 p1[i]=p1[i]*1.0e9;
 cout<<"Get parameter [0] is "<<p0[i]<<endl;
 cout<<"Get parameter [1] is "<<p1[i]<<endl;

  fx[i]->Draw("same"); 
  cxp[i]->cd(3);
 
 //----- Corrected for X' -------------//  
 
 hxp_c[i] = new TH2F(Form("hxp_c[%d]",i),Form("pad #%d RF -s2 time(s2 X' corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.05,0.0);

 T->Project(Form("hxp_c[%d]",i),Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0-R.tr.th[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i]);


 hxp_c[i]->Draw();
 cxp[i]->cd(4);
 

 htof_xpc[i] = new TH1F(Form("htof_xpc[%d]",i),Form("pad #%d RF -s2 time(s2 X' corrected)",i),300,0.2e-6,0.25e-6);

 T->Project(Form("htof_xpc[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -R.tr.th[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i]&& cut_tof[i]);

 htof[i]->Draw();
 htof_xpc[i]->SetLineColor(1);
 htof_xpc[i]->Draw("same");




 /*
 //============= Y correction ============================//



 //========== Y' Correction ====================//

cyp[i] = new TCanvas(Form("cyp[%d]",i),Form("cyp[%d]",i));
   cyp[i]->Divide(2,2);


//-------- TOF vs X' hist --------------------------//
 hyp[i] = new TH2F(Form("hyp[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.05,0.0);

T->Project(Form("hyp[%d]",i),Form("(R.tr.ph[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",pmtl+i,pmtr+i),cut[i]);
cyp[i]->cd(1);
 hyp[i]->Draw();


//---------TProfile  Correction for Y' ----------------------//
 cyp[i]->cd(2);

 tpyp[i] = new TProfile(Form("tpyp[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;s2 Y' position  ; ",i),200,-0.05,0.00,0.2*1.0e-6,0.25*1.0e-6);

 T->Project(Form("tpyp[%d]",i),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(R.tr.ph[0])",pmtl+i,pmtr+i),cut[i]);
 tpyp[i]->Draw();
 //------ Fit TProfile Correction for Y' -------------------//

 fyp[i]= new TF1(Form("fyp[%d]",i),"[0]+[1]*x",0.2206e-6,0.224e-6);
 tpyp[i]->Fit(Form("fyp[%d]",i),"","");
 p0[i]=fyp[i]->GetParameter(0);
 p1[i]=fyp[i]->GetParameter(1);

 p1[i]=p1[i]*1.0e9;
 cout<<"Get parameter [0] is "<<p0[i]<<endl;
 cout<<"Get parameter [1] is "<<p1[i]<<endl;

  fy[i]->Draw("same"); 
  cyp[i]->cd(3);
 
 //----- Corrected for Y' -------------//  
 
 hyp_c[i] = new TH2F(Form("hyp_c[%d]",i),Form("pad #%d RF -s2 time(s2 Y' corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,-0.05,0.0);

 T->Project(Form("hyp_c[%d]",i),Form("(R.tr.ph[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0-R.tr.ph[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i]);


 hyp_c[i]->Draw();
 cyp[i]->cd(4);
 

 htof_ypc[i] = new TH1F(Form("htof_ypc[%d]",i),Form("pad #%d RF -s2 time(s2 Y' corrected)",i),300,0.2e-6,0.25e-6);

 T->Project(Form("htof_ypc[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -R.tr.ph[0]*%lf*1.0e-9)",pmtl+i,pmtr+i,p1[i]),cut[i]&& cut_tof[i]);

 htof[i]->Draw();
 htof_ypc[i]->SetLineColor(1);
 htof_ypc[i]->Draw("same");


 */





   }
 }


 /*  
  //// X vs RF -s2 (S2 is corrected for X)
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Draw(Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.25,0.2235e-6,-0.02);
 l1->SetLineColor(kRed);
  l1->Draw();


 

  /// X' vs RF -s2 now s2 is already corrected for X

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.05,0);

 T->Draw(Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.05,0.2235e-6,0);
 l1->SetLineColor(kRed);
  l1->Draw();
  
  /// Correction for X'

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.05,0,300,0.2218e-6,0.2252e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  ) :(R.tr.th[0])>>h8",pmt+15,pmt+47),cut_8,"colz");
 // For TProfile
TProfile *ht8 = new TProfile("ht8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,-0.05,0,0.2218e-6,0.2252e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  ) :(R.tr.th[0])>>ht8",pmt+15,pmt+47),cut_8,"same");


TF1 *ft = new TF1("ft","[0]+[1]*x",-0.037,-0.015);
  ht8->Fit("ft","R+");
  
  // Checking for X' vs RF -s2

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.05,0);

 T->Draw(Form("(R.tr.th[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0]  + 4.27527e-08*R.tr.th[0])>>h8",pmt+15,pmt+47),cut_8,"colz");
 //   p0           2.22662e-07 
  TLine *l1 = new TLine(0.2228e-6,-0.05,0.2228e-6,0);
 l1->SetLineColor(kRed);
  l1->Draw();

  
  // just check the X vs RF -s2

 TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.25,-0.02);

 T->Draw(Form("(R.tr.x[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.23500e-07 
  TLine *l1 = new TLine(0.2235e-6,-0.25,0.2235e-6,-0.02);
 l1->SetLineColor(kRed);
  l1->Draw();


  // ADCL vs RF -s2 

 TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,120,350);

 T->Draw(Form("(R.s2.la_c[7]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");
 
  TLine *l1 = new TLine(0.2235e-6,120,0.2235e-6,350);
 l1->SetLineColor(kRed);
  l1->Draw();
  // ADCL vs RF -S2

TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,120,350,300,0.221e-6,0.2235e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] ):(R.s2.la_c[7]) >>h8",pmt+15,pmt+47),cut_8,"colz");
 
 // For tProfile
TProfile *ht8 = new TProfile("ht8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",200,120,350,0.221e-6,0.2235e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] ):(R.s2.la_c[7]) >>ht8",pmt+15,pmt+47),cut_8,"same");
 
TF1 *ft = new TF1("ft","[0]+[1]*x",185,250);
  ht8->Fit("ft","R+");
  
  // checking ADCL vs RF -s2 after ADC correction


TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,120,350);

 T->Draw(Form("(R.s2.ra_c[7]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] -2.37582e-13*R.s2.la_c[7]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 //  p0           2.22499e-07 
  TLine *l1 = new TLine(0.2235e-6,120,0.2235e-6,350);
 l1->SetLineColor(kRed);
  l1->Draw();


  /// R.tr.y[0]  vs RF-s2
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH2F *h8 = new TH2F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6,200,-0.04,0.05);

 T->Draw(Form("(R.tr.ph[0]):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] -2.37582e-13*R.s2.la_c[7]  )>>h8",pmt+15,pmt+47),cut_8,"colz");
 // This is final Histogram counts vs RF -s2 
TCut cut_8 =("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==7    &&  R.s2.nthit==1");

 TCanvas *c8 = new TCanvas("c8","c8", 600,600);
 
 TH1F *h8 = new TH1F("h8","pad #7 RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",300,0.2e-6,0.25e-6);

 T->Draw(Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -8.23538e-09*R.tr.x[0] + 4.27527e-08*R.tr.th[0] )>>h8",pmt+15,pmt+47),cut_8,"colz");

 //+2.37582e-13*R.s2.la_c[7]
  

 */
}




//========== Cut Condition ============//

/*
Tcut cut_tof(int i){

   int chmax=26;
  TCut cut_tof;
  double min[chmax],max[chmax];





  // ====== peak position ============= //
  for(int j=0;j<chmax;j++){
    min[j]=0;
    max[j]=0;
  } 


  min[7]=0.225;
  max[7]=0.227;


    cut_tof=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>%lf*1.0e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<%lf*1e-6",pmtl+i,pmtr+i,min[i],pmtl+i,pmtr+i,max[i]);





 cut_tof[4]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>0.225e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<0.225e-6",pmtl+4,pmtr+4,pmtl+4,pmtr+4);
  }{}
 cut_tof[5]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>0.225e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<0.225e-6",pmtl+5,pmtr+5,pmtl+5,pmtr+5);

  cut_tof[7]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>0.225e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<0.227e-6 ",pmtl+7,pmtr+7,pmtl+7,pmtr+7);

 cut_tof[8]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>0.2205e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<0.222e-6 ",pmtl+8,pmtr+8,pmtl+8,pmtr+8);
   

 return cut_tof; 

}

*/
