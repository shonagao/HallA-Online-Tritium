


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

void ehodo_TOF_correct1330(){
  gStyle->SetOptFit(111);
  
  // ===================================== //
  // ====== Open a ROOT file ============= //
 // ===================================== //
  TFile* f1 = new TFile("hkstest_01330.root"); // position was displaed from center
  TTree* t1 = (TTree*)f1->Get("T");
  
  // ===================================== //
  // ====== Difine cut conditions ======== //
  // ===================================== //
  TCut cut_tdc00 = "D.TH00[0]>600 && D.TH00[0]<1200 && Ndata.D.TH00==1";
  TCut cut_tdc01 = "D.TH01[0]>300 && D.TH01[0]<800 && Ndata.D.TH01==1";
  TCut cut_tdc02 = "D.TH02[0]>400 && D.TH02[0]<1200 && Ndata.D.TH02==1";
  TCut cut_tdc03 = "D.TH03[0]>700 && D.TH03[0]<800 && Ndata.D.TH03==1";

  TCut cut_qdc00 = "D.A00>900 && D.A00<1400";
  TCut cut_qdc01 = "D.A01>800 && D.A01<1000";
  TCut cut_qdc02 = "D.A02>1200 && D.A02<1700";
  TCut cut_qdc03 = "D.A03>1000 && D.A03<1500";
  TCut cut_other ="-100<(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2 && 200>(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2";
  TCut cutall = "cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 && cut_other";
  

  // ================================== //
  // ====== Create histograms ========= //
  // ================================== //

  /*

  ////////D.A hist NOT cut condition hist///////////////////////
  TCanvas *cA = new TCanvas("cA","cA");
  //////////hist water chelencov //////////////////////
 int ybin = 1200;
  double ymin = -200.0, ymax = 1000.0;
  TH1F* histA04 = new TH1F("histA04","",ybin,ymin,ymax);
  t1->Project("histA04","D.A04");
  
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  
 

  histA04->GetXaxis()->SetRangeUser(600,800);
  histA04->GetXaxis()->SetTitle("D.A04 [ch]");
  histA04->GetYaxis()->SetTitle("Counts / 50 ps");
  histA04->Draw();
 
TCanvas *cB = new TCanvas("cB","cB");
TH2F* histB = new TH2F("histB","",100,-500.0,1500.0,100,-500.0,1500);
  t1->Project("histB",
	      "D.TH00:D.A00"
	      ,
	      //cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 &&
	       cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      //&& cut_other
);


  histB->GetXaxis()->SetTitle("ADC [ch]");
  histB->GetYaxis()->SetTitle("TDC [ch]");
  //hist00->GetYaxis()->SetRangeUser(-20,6);
  histB->Draw("colz");
 

*/

  ////////TOF hist before correction//////////////
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


  /*
  //===================ADC:TDC=====================================//
  TCanvas *c3 =new TCanvas("c3","c3");
  c3->Divide(4,2);
  c3->cd(1);
  TH2F* hist20 =new TH2F("hist20","",100,600.0,1500.0,xbin*10,xmin*10,xmax*15);
    t1->Project("hist20","D.TH00:D.A00"
		//	,cutall

		);

    hist20->Draw("colz");
  c3->cd(2);
  TH2F* hist21 =new TH2F("hist21","",100,600.0,1500.0,xbin*10,xmin*10,xmax*15);
    t1->Project("hist21","D.TH01:D.A00"
		//	,cutall

		);

 hist21->Draw("colz");
 c3->cd(3);
 TH2F* hist22 =new TH2F("hist22","",100,600.0,1500.0,xbin*10,xmin*10,xmax*15);
    t1->Project("hist22","D.TH02:D.A00"
		//	,cutall

		);

 hist22->Draw("colz");
 c3->cd(4);
 TH2F* hist23 =new TH2F("hist23","",100,600.0,1500.0,xbin*10,xmin*10,xmax*15);
    t1->Project("hist23","D.TH03:D.A00"
		//	,cutall

		);
 hist23->Draw("colz");
*/
  //==========TOF:ADC==============================//

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(4,2);

  c0->cd(1);

TH2F* hist00 = new TH2F("hist00","",100,600.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00"
	      ,
	        cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03   && 
           cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	   &&    cut_other);


  hist00->GetXaxis()->SetTitle("ADC0 [ch]");
  hist00->GetYaxis()->SetTitle("TOF [ns]");
  hist00->GetXaxis()->SetRangeUser(800,1400);
  hist00->GetYaxis()->SetRangeUser(2,4);
  hist00->Draw("colz");
  c0->cd(2);

TH2F* hist01 = new TH2F("hist01","",100,600.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 &&
	      cut_qdc00 && cut_qdc01
	      && cut_qdc02 && cut_qdc03 
	      &&
 cut_other);
  hist01->GetXaxis()->SetTitle("ADC1 [ch]");
  hist01->GetYaxis()->SetTitle("TOF [ns]");
  hist01->GetXaxis()->SetRangeUser(800,1100);
  hist01->GetYaxis()->SetRangeUser(2,4);
  hist01->Draw("colz");
  c0->cd(3);

TH2F* hist02 = new TH2F("hist02","",100,600.0,1700.0,xbin,xmin,xmax);
  t1->Project("hist02",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A02",
	           cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 &&
	      cut_qdc00 && cut_qdc01 && cut_qdc02 &&
	      cut_qdc03 	      && cut_other);
  hist02->GetXaxis()->SetTitle("ADC2 [ch]");
  hist02->GetYaxis()->SetTitle("TOF [ns]");

  hist02->GetYaxis()->SetRangeUser(2,4);
 hist02->Draw("colz");
  c0->cd(4);
hist02->GetXaxis()->SetRangeUser(1200,1700);
TH2F* hist03 = new TH2F("hist03","",100,600.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist03",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A03",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 &&
	       cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	     && cut_other);
  hist03->GetXaxis()->SetTitle("ADC3 [ch]");
  hist03->GetYaxis()->SetTitle("TOF [ns]");
hist03->GetXaxis()->SetRangeUser(1000,1500);
  hist03->GetYaxis()->SetRangeUser(2,4);
 hist03->Draw("colz");
 c0->cd(5);

  // ============================== //
  // ===== TPlofile ================ //
  // ============================= //
 //  c0->cd(1);
TProfile *pf00 = new TProfile("pf00","",50,700,1700);
  t1->Project("pf00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  //pf00->Draw("same");
  // c0->cd(2);

TProfile *pf01 = new TProfile("pf01","",50,700,1700);
  t1->Project("pf01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  // pf01->Draw("same");
  //c0->cd(3);

TProfile *pf02 = new TProfile("pf02","",50,700,1700);
  t1->Project("pf02",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A02",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  //pf02->Draw("same");
  //c0->cd(4);

TProfile *pf03 = new TProfile("pf03","",50,700,1700);
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


  TF1* ffunc00 = new TF1("ffunc00","[0]/sqrt(x)+[1]",500,2000);
 ffunc00->SetNpx(1000);
  pf00->Fit("ffunc00","","",900,1200);

    Double_t p00 = ffunc00->GetParameter(0);
    Double_t q00 = ffunc00->GetParameter(1);

    //ffunc00->Draw("same");
    // c0->cd(2);

  TF1* ffunc01 = new TF1("ffunc01","[0]/sqrt(x)+[1]",500,2000);
 ffunc01->SetNpx(1000);
  pf01->Fit("ffunc01","","",850,1000);

    Double_t p01 = ffunc01->GetParameter(0);
    Double_t q01 = ffunc01->GetParameter(1);

    // ffunc01->Draw("same");
    //c0->cd(3);

  TF1* ffunc02 = new TF1("ffunc02","[0]/sqrt(x)+[1]",500,1000);
 ffunc02->SetNpx(1000);
  pf02->Fit("ffunc02","","",1400,1600);

    Double_t p02 = ffunc02->GetParameter(0);
    Double_t q02 = ffunc02->GetParameter(1);

    //ffunc02->Draw("same");
    // c0->cd(4);

  TF1* ffunc03 = new TF1("ffunc03","[0]/sqrt(x)+[1]",500,2000);
 ffunc03->SetNpx(1000);
  pf03->Fit("ffunc03","","",1100,1300);

    Double_t p03 = ffunc03->GetParameter(0);
    Double_t q03 = ffunc03->GetParameter(1);

    // ==================================== //
    // ====== Correction ================= //
    // ================================== //


    //ffunc03->Draw("same");

    //c0->cd(5);
  
TH2F* hist10 = new TH2F("hist10","",100,600.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist10",
   Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A00",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
 
  hist10->GetXaxis()->SetRangeUser(800,1400);
   hist10->GetYaxis()->SetRangeUser(-10,-8);
 hist10->GetXaxis()->SetTitle("ADC0 [ch]");
  hist10->GetYaxis()->SetTitle("TOF [ns]");
 hist10->Draw("colz");
  c0->cd(6);

TH2F* hist11 = new TH2F("hist11","",100,600.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist11",
	       Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A01",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
 hist11->GetXaxis()->SetTitle("ADC1 [ch]");
  hist11->GetYaxis()->SetTitle("TOF [ns]");
hist11->GetXaxis()->SetRangeUser(800,1100);
 hist11->GetYaxis()->SetRangeUser(-10,-8);
 hist11->Draw("colz");
  c0->cd(7);

TH2F* hist12 = new TH2F("hist12","",100,600.0,1700.0,xbin,xmin,xmax);
  t1->Project("hist12",
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A02",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  hist12->GetXaxis()->SetTitle("ADC2 [ch]");
  hist12->GetYaxis()->SetTitle("TOF [ns]");
  hist12->GetXaxis()->SetRangeUser(1200,1700);
  hist12->GetYaxis()->SetRangeUser(-10,-8);
  hist12->Draw("colz");
  c0->cd(8);


TH2F* hist13 = new TH2F("hist13","",100,500.1000,1500.0,xbin,xmin,xmax);
  t1->Project("hist13",
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)-(%lf/sqrt(D.A01)+%lf)-(%lf/sqrt(D.A02)+%lf)-(%lf/sqrt(D.A03)+%lf):D.A03",p00,q00,p01,q01,p02,q02,p03,q03),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
    
 hist13->GetXaxis()->SetTitle("ADC3 [ch]");
 hist13->GetYaxis()->SetTitle("TOF [ns]");
 hist13->GetXaxis()->SetRangeUser(1000,1500);
 hist13->GetYaxis()->SetRangeUser(-10,-8);
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
  hist1->Fit("ffunc1","","",0.0,4.0);
  hist1->GetXaxis()->SetRangeUser(0,6);
  hist1->SetTitleSize(0.05,"x");
  hist1->SetTitleSize(0.05,"y");
  hist1->GetXaxis()->CenterTitle();
  hist1->GetYaxis()->CenterTitle();
   
  // ================================== //
  // ====== Draw the histograms ======= //
  // ================================== //
 
    hist1->Draw();





  //-----------After Correction-----------//

    // TCanvas *c2 = new TCanvas("c2","c2");
   c1->cd(2);
  
 TF1* ffunc5 = new TF1("ffunc5","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-10,0.0);
 ffunc5->SetNpx(1000);
 ffunc5->SetParameters(2000.5,-9.5,2);
 hist5->Fit("ffunc5","","",-12.5,-8);
 hist5->GetXaxis()->SetRangeUser(-12,-6);
  hist5->SetTitleSize(0.05,"x");
  hist5->SetTitleSize(0.05,"y");
  hist5->GetXaxis()->CenterTitle();
  hist5->GetYaxis()->CenterTitle();
// ================================== //
  // ====== Draw the histograms ======= //
  // ================================== //
 
   hist5->Draw();

TCanvas *c2 = new TCanvas("c2","c2");
 c2->Divide(2,1); 
c2->cd(1); 
hist00->Draw("colz");
  hist00->GetXaxis()->SetTitle("ADC0 [ch]");
  hist00->GetYaxis()->SetTitle("TOF [ns]");
  hist00->GetXaxis()->SetRangeUser(800,1400);
  hist00->GetYaxis()->SetRangeUser(2,4);
  hist00->SetTitleSize(0.05,"x");
  hist00->SetTitleSize(0.05,"y");
  hist00->GetXaxis()->CenterTitle();
  hist00->GetYaxis()->CenterTitle();
hist00->Draw("colz");


  pf00->Draw("same");
  ffunc00->Draw("same");
 c2->cd(2);
 
hist10->Draw("colz");
  hist10->GetXaxis()->SetRangeUser(800,1400);
   hist10->GetYaxis()->SetRangeUser(-10,-8);
 hist10->GetXaxis()->SetTitle("ADC0 [ch]");
  
hist10->SetTitleSize(0.05,"x");
  hist10->SetTitleSize(0.05,"y");
hist10->GetYaxis()->SetTitle("TOF [ns]");
  hist10->GetXaxis()->CenterTitle();
  hist10->GetYaxis()->CenterTitle(); 
hist10->Draw("colz");


}













