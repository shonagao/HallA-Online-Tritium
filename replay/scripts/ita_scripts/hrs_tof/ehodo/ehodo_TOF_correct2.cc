
/*
  "gtest.C"
  Quick check of TOF resolution (EHODO)
  HV = -1700 V for all channels
  
  Modified from test.C coded by K.Itabashi
  Nov 14, 2017


  "ehodo_TOF_correction.cc"
  Modefied gtest.C coded by K.Itabashi
  Nov 16, 2017

  "ehodo_TOF_correction2.cc"
  Modefied ehodo_TOF_correctioncc  coded by K.Itabashi
  Nov 20,2017
  run1321

 

*/

void ehodo_TOF_correct2(){
  gStyle->SetOptFit(111);
  
  // ===================================== //
  // ====== Open a ROOT file ============= //
 // ===================================== //
  TFile* f1 = new TFile("hkstest_01321.root"); // position was displaed from center
  TTree* t1 = (TTree*)f1->Get("T");
  
  // ===================================== //
  // ====== Difine cut conditions ======== //
  // ===================================== //
  TCut cut_tdc00 = "D.TH00[0]>620 && D.TH00[0]<650 && Ndata.D.TH00==1";
  TCut cut_tdc01 = "D.TH01[0]>350 && D.TH01[0]<500 && Ndata.D.TH01==1";
  TCut cut_tdc02 = "D.TH02[0]>400 && D.TH02[0]<600 && Ndata.D.TH02==1";
  TCut cut_tdc03 = "D.TH03>350 && D.TH03<550 && Ndata.D.TH03==1";
  TCut cut_qdc00 = "D.A00>500 && D.A00<1000";
  TCut cut_qdc01 = "D.A01>500 && D.A01<1000";
  TCut cut_qdc02 = "D.A02>600 && D.A02<1200";
  TCut cut_qdc03 = "D.A03>500 && D.A03<1000";
  TCut cut_other ="-50<(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2 && 100>(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2";
  TCut cutall = "cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 && cut_other";
  

  // ================================== //
  // ====== Create histograms ========= //
  // ================================== //
  int xbin = 1200;
  double xmin = -20.0, xmax = 100.0;
  TH1F* hist1 = new TH1F("hist1","",xbin,xmin,xmax);
  t1->Project("hist1",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  hist1->GetXaxis()->SetRangeUser(-6,6);
  hist1->GetXaxis()->SetTitle("TOF [ns]");
  hist1->GetYaxis()->SetTitle("Counts / 50 ps");

 


  //==========TOF:ADC==============================//

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(4,2);

  c0->cd(1);

TH2F* hist00 = new TH2F("hist00","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  hist00->GetXaxis()->SetTitle("ADC0 [ch]");
  hist00->GetYaxis()->SetTitle("TOF [ns]");
hist00->GetYaxis()->SetRangeUser(-20,6);
  hist00->Draw("colz");
  c0->cd(2);

TH2F* hist01 = new TH2F("hist01","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  hist01->GetXaxis()->SetTitle("ADC1 [ch]");
  hist01->GetYaxis()->SetTitle("TOF [ns]");
hist01->GetYaxis()->SetRangeUser(-20,6);
  hist01->Draw("colz");
  c0->cd(3);

TH2F* hist02 = new TH2F("hist02","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist02",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A02",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  hist02->GetXaxis()->SetTitle("ADC2 [ch]");
  hist02->GetYaxis()->SetTitle("TOF [ns]");
hist02->GetYaxis()->SetRangeUser(-20,6);
 hist02->Draw("colz");
  c0->cd(4);

TH2F* hist03 = new TH2F("hist03","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist03",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A03",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  hist03->GetXaxis()->SetTitle("ADC3 [ch]");
  hist03->GetYaxis()->SetTitle("TOF [ns]");
hist03->GetYaxis()->SetRangeUser(-20,6);
 hist03->Draw("colz");
 c0->cd(5);
  // ============================== //
  // ===== TPlofile ================ //
  // ============================= //
  //c0->cd(1);
TProfile *pf00 = new TProfile("pf00","",50,500,1000);
  t1->Project("pf00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  //pf00->Draw("same");
  // c0->cd(2);

TProfile *pf01 = new TProfile("pf01","",50,500,1000);
  t1->Project("pf01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  // pf01->Draw("same");
  //c0->cd(3);

TProfile *pf02 = new TProfile("pf02","",50,500,1000);
  t1->Project("pf02",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A02",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  //pf02->Draw("same");
  //c0->cd(4);

TProfile *pf03 = new TProfile("pf03","",50,500,1000);
  t1->Project("pf03",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A03",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  //pf03->Draw("same");
  //c0->cd(1);

  // =========================================== //
  // ========== Get Function Parameter ======== //
  // ========================================= //


  TF1* ffunc00 = new TF1("ffunc00","[0]/sqrt(x)+[1]",500,1000);
 ffunc00->SetNpx(1000);
  pf00->Fit("ffunc00","","",600,800);

    Double_t p00 = ffunc00->GetParameter(0);
    Double_t q00 = ffunc00->GetParameter(1);

    //ffunc00->Draw("same");
    // c0->cd(2);

  TF1* ffunc01 = new TF1("ffunc01","[0]/sqrt(x)+[1]",500,1000);
 ffunc01->SetNpx(1000);
  pf01->Fit("ffunc01","","",600,750);

    Double_t p01 = ffunc01->GetParameter(0);
    Double_t q01 = ffunc01->GetParameter(1);

    // ffunc01->Draw("same");
    //c0->cd(3);

  TF1* ffunc02 = new TF1("ffunc02","[0]/sqrt(x)+[1]",500,1000);
 ffunc02->SetNpx(1000);
  pf02->Fit("ffunc02","","",700,850);

    Double_t p02 = ffunc02->GetParameter(0);
    Double_t q02 = ffunc02->GetParameter(1);

    //ffunc02->Draw("same");
    // c0->cd(4);

  TF1* ffunc03 = new TF1("ffunc03","[0]/sqrt(x)+[1]",500,1000);
 ffunc03->SetNpx(1000);
  pf03->Fit("ffunc03","","",600,700);

    Double_t p03 = ffunc03->GetParameter(0);
    Double_t q03 = ffunc03->GetParameter(1);

    // ==================================== //
    // ====== Correction ================= //
    // ================================== //


    //ffunc03->Draw("same");

    //c0->cd(5);
  
TH2F* hist10 = new TH2F("hist10","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist10",
   Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A00",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
   hist10->GetYaxis()->SetRangeUser(-20,6);
 hist10->GetXaxis()->SetTitle("ADC0 [ch]");
  hist10->GetYaxis()->SetTitle("TOF [ns]");
 hist10->Draw("colz");
  c0->cd(6);

TH2F* hist11 = new TH2F("hist11","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist11",
	       Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A01",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
 hist11->GetXaxis()->SetTitle("ADC1 [ch]");
  hist11->GetYaxis()->SetTitle("TOF [ns]");
 hist11->GetYaxis()->SetRangeUser(-20,6);
 hist11->Draw("colz");
  c0->cd(7);

TH2F* hist12 = new TH2F("hist12","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist12",
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A02",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
 hist12->GetXaxis()->SetTitle("ADC2 [ch]");
  hist12->GetYaxis()->SetTitle("TOF [ns]");
   hist12->GetYaxis()->SetRangeUser(-20,6);
hist12->Draw("colz");
  c0->cd(8);


TH2F* hist13 = new TH2F("hist13","",100,500.0,1000.0,xbin,xmin,xmax);
  t1->Project("hist13",
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A03",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
    
 hist13->GetXaxis()->SetTitle("ADC3 [ch]");
  hist13->GetYaxis()->SetTitle("TOF [ns]");
 hist13->GetYaxis()->SetRangeUser(-20,6);
hist13->Draw("colz");


//  TCanvas *cx = new TCanvas("cx","cx");

TH1F* hist5 = new TH1F("hist5","",xbin,xmin,xmax);
  t1->Project("hist5",
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf)",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  hist5->GetXaxis()->SetRangeUser(-200.0,100.0);
  hist5->GetXaxis()->SetTitle("TOF [ns]");
  hist5->GetYaxis()->SetTitle("Counts / 50 ps");

  //hist5->Draw("colz");


  // ==================================== //
  // ===== Fitting to the spectrom ====== //
  // ==================================== //
   
  //------------Before Correction---------//

    TCanvas *c1 = new TCanvas("c1","c1");
   c1->cd();
   c1->Divide(2,1);
   c1->cd(1);

  TF1* ffunc1 = new TF1("ffunc1","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-5.0,5.0);
  ffunc1->SetNpx(1000);
  ffunc1->SetParameters(1000,3,2);
  hist1->Fit("ffunc1","","",0.0,10.0);
  //hist1->GetYaxis()->SetRangeUser(0,6);

   
  // ================================== //
  // ====== Draw the histograms ======= //
  // ================================== //
 
    hist1->Draw();





  //-----------After Correction-----------//

    // TCanvas *c2 = new TCanvas("c2","c2");
   c1->cd(2);
  
 TF1* ffunc5 = new TF1("ffunc5","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-10,0.0);
 ffunc5->SetNpx(1000);
 ffunc5->SetParameters(1000,-9.0,0.1);
 hist5->Fit("ffunc5","","",-10.0,0.0);
hist5->GetXaxis()->SetRangeUser(-10,0);
// ================================== //
  // ====== Draw the histograms ======= //
  // ================================== //
 
   hist5->Draw();

}

