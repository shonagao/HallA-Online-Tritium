

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

void ehodo_TOF_correctnew1330(){
  gStyle->SetOptFit(111);
  
  // ===================================== //
  // ====== Open a ROOT file ============= //
 // ===================================== //
  TFile* f1 = new TFile("hkstest_01330.root"); // position was displaed from center
  TTree* t1 = (TTree*)f1->Get("T");
  
  // ===================================== //
  // ====== Difine cut conditions ======== //
  // ===================================== //
  TCut cut_tdc00 = "D.TH00[0]>300 && D.TH00[0]<2000 && Ndata.D.TH00==1";
  TCut cut_tdc01 = "D.TH01[0]>300 && D.TH01[0]<2000 && Ndata.D.TH01==1";
  TCut cut_tdc02 = "D.TH02[0]>300 && D.TH02[0]<2000 && Ndata.D.TH02==1";
  TCut cut_tdc03 = "D.TH03[0]>300 && D.TH03[0]<2000 && Ndata.D.TH03==1";

  TCut cut_qdc00 = "D.A00>500 && D.A00<2200";
  TCut cut_qdc01 = "D.A01>500 && D.A01<2200";
  TCut cut_qdc02 = "D.A02>500 && D.A02<2200";
  TCut cut_qdc03 = "D.A03>500 && D.A03<2200";
  TCut cut_other ="-100<(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2 && 200>(D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2";
  TCut cutall = "cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 && cut_other";
  

  // ================================== //
  // ====== Create histograms ========= //
  // ================================== //

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

  //===================ADC:TDC=====================================//

  int ybin =200;
  double ymin =30, ymax =50;
  TCanvas *c3 =new TCanvas("c3","c3");
  c3->Divide(4,2);
  c3->cd(1);
  TH2F* hist20 =new TH2F("hist20","",100,600.0,1500.0,ybin,ymin,ymax);
    t1->Project("hist20","D.TH00*0.05:D.A00"
			, cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other

		);
  hist20->GetXaxis()->SetTitle("ADC0 [ch]");
  hist20->GetYaxis()->SetTitle("TDC0[/50 ps]");
  hist20->GetXaxis()->SetRangeUser(500,1200);
  //hist20->GetYaxis()->SetRangeUser(45,47);
 
    hist20->Draw("colz");
  c3->cd(2);
  TH2F* hist21 =new TH2F("hist21","",100,600.0,1500.0,ybin,ymin,ymax);
    t1->Project("hist21","D.TH01*0.05:D.A01"
    			, cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other


		);
  hist21->GetXaxis()->SetTitle("ADC1 [ch]");
  hist21->GetYaxis()->SetTitle("TDC1[ch]");
  hist21->GetXaxis()->SetRangeUser(500,1200);
  //  hist21->GetYaxis()->SetRangeUser(35,37);
 hist21->Draw("colz");
 c3->cd(3);
 TH2F* hist22 =new TH2F("hist22","",100,600.0,1700.0,ybin,ymin,ymax);
    t1->Project("hist22","D.TH02*0.05:D.A02"
 			, cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other


		);
  hist22->GetXaxis()->SetTitle("ADC2 [ch]");
  hist22->GetYaxis()->SetTitle("TDC2[ch]");
  hist22->GetXaxis()->SetRangeUser(500,1600);
  // hist22->GetYaxis()->SetRangeUser(38,40);
 hist22->Draw("colz");
 c3->cd(4);
 TH2F* hist23 =new TH2F("hist23","",100,600.0,1700.0,ybin,ymin,ymax);
    t1->Project("hist23","D.TH03*0.05:D.A03"			, cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other


		);
  hist23->GetXaxis()->SetTitle("ADC3 [ch]");
  hist23->GetYaxis()->SetTitle("TDC3[ch]");
  hist23->GetXaxis()->SetRangeUser(500,1300);
  //hist23->GetYaxis()->SetRangeUser(36,38);
 hist23->Draw("colz");

  //==========TOF:ADC==============================//

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(4,2);

  c0->cd(1);

TH2F* hist00 = new TH2F("hist00","",100,600.0,1700.0,xbin,xmin,xmax);
  t1->Project("hist00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00"
	      ,
	        cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03   && 
           cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	   &&    cut_other);


  hist00->GetXaxis()->SetTitle("ADC0 [ch]");
  hist00->GetYaxis()->SetTitle("TOF [ns]");
  hist00->GetXaxis()->SetRangeUser(800,1400);
  // hist00->GetYaxis()->SetRangeUser(2,4);
  hist00->Draw("colz");
  c0->cd(2);

TH2F* hist01 = new TH2F("hist01","",100,500.0,1500.0,xbin,xmin,xmax);
  t1->Project("hist01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 &&
	      cut_qdc00 && cut_qdc01
	      && cut_qdc02 && cut_qdc03 
	      &&
 cut_other);
  hist01->GetXaxis()->SetTitle("ADC1 [ch]");
  hist01->GetYaxis()->SetTitle("TOF [ns]");
  hist01->GetXaxis()->SetRangeUser(500,1100);
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



 //===================================================//
 //======TDC:QDC TPlofire=============================//
 //===================================================//


   c3->cd(1);
TProfile *pf20 = new TProfile("pf20","",50,900,1200);
  t1->Project("pf20",
	      "D.TH00*0.05:D.A00",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
   pf20->Draw("same");
   c3->cd(2);

TProfile *pf21 = new TProfile("pf21","",50,700,1700);
  t1->Project("pf21",
	      "D.TH01*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  pf21->Draw("same");
  c3->cd(3);

TProfile *pf22 = new TProfile("pf22","",50,700,1700);
  t1->Project("pf22",
	      "D.TH02*0.05:D.A02",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

   pf22->Draw("same");
  c3->cd(4);

TProfile *pf23 = new TProfile("pf23","",50,700,1700);
  t1->Project("pf23",
	      "D.TH03*0.05:D.A03",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  pf23->Draw("same");
  c0->cd(1);















  //=====================================================//
  //=================TOF:QDC TProfile============================//
  //===================================================//



   c0->cd(1);
TProfile *pf00 = new TProfile("pf00","",50,700,1700);
  t1->Project("pf00",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A00",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
  pf00->Draw("same");
   c0->cd(2);

TProfile *pf01 = new TProfile("pf01","",50,700,1700);
  t1->Project("pf01",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A01",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
  
   pf01->Draw("same");
  c0->cd(3);

TProfile *pf02 = new TProfile("pf02","",50,700,1700);
  t1->Project("pf02",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A02",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  pf02->Draw("same");
  c0->cd(4);

TProfile *pf03 = new TProfile("pf03","",50,700,1700);
  t1->Project("pf03",
	      "((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05:D.A03",
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);

  pf03->Draw("same");
  c0->cd(1);

  // =========================================== //
  // ========== Get Function Parameter ======== //
  // ========================================= //

  //=================================================================//
  //===================TDC:QDC Fitting=============================//
  //===============================================================//

  c3->cd(1);
 TF1* ffunc20 = new TF1("ffunc20","[0]/sqrt(x)+[1]",900,1300);
 ffunc20->SetNpx(1000);
  pf20->Fit("ffunc20","","",900,1200);

    Double_t p20 = ffunc20->GetParameter(0);
    Double_t q20 = ffunc20->GetParameter(1);

    ffunc20->Draw("same");
     c3->cd(2);

  TF1* ffunc21 = new TF1("ffunc21","[0]/sqrt(x)+[1]",500,2000);
 ffunc21->SetNpx(1000);
  pf21->Fit("ffunc21","","",850,1000);

    Double_t p21 = ffunc21->GetParameter(0);
    Double_t q21 = ffunc21->GetParameter(1);

    ffunc21->Draw("same");
    c3->cd(3);

  TF1* ffunc22 = new TF1("ffunc22","[0]/sqrt(x)+[1]",500,1000);
 ffunc22->SetNpx(1000);
  pf22->Fit("ffunc22","","",1300,1600);

    Double_t p22 = ffunc22->GetParameter(0);
    Double_t q22 = ffunc22->GetParameter(1);

      ffunc22->Draw("same");
       c3->cd(4);

  TF1* ffunc23 = new TF1("ffunc23","[0]/sqrt(x)+[1]",500,2000);
 ffunc23->SetNpx(1000);
  pf23->Fit("ffunc23","","",1000,1300);

    Double_t p23 = ffunc23->GetParameter(0);
    Double_t q23 = ffunc23->GetParameter(1);

    ffunc23->Draw("same");
    c0->cd(1);

  //==============================================//
  //=======TOF:QDC fitting=======================//
  //=============================================//

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


    //==========TDC:QDC correction================//




    c3->cd(5);
TH2F* hist30 = new TH2F("hist30","",100,600.0,1500.0,ybin*10,0,100);
  t1->Project("hist30",
   Form("(D.TH00-(%lf/sqrt(D.A00)+%lf))*0.05:D.A00",p20,q20),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
 
  hist30->GetXaxis()->SetRangeUser(890,1200);
  hist30->GetYaxis()->SetRangeUser(40,50);
  hist30->GetXaxis()->SetTitle("ADC0 [ch]");
  hist30->GetYaxis()->SetTitle("TDC0 [ch/50 ps]");
  hist30->Draw("colz");
  c3->cd(6);




TH2F* hist31 = new TH2F("hist31","",100,600.0,1500.0,ybin*10,0,100);
  t1->Project("hist31",
   Form("(D.TH01-(%lf/sqrt(D.A01)+%lf))*0.05:D.A01",p21,q21),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
 
  hist31->GetXaxis()->SetRangeUser(900,1200);
  hist31->GetYaxis()->SetRangeUser(30,40);
  hist31->GetXaxis()->SetTitle("ADC1 [ch]");
  hist31->GetYaxis()->SetTitle("TDC1 [ch/50 ps]");
  hist31->Draw("colz");


    c3->cd(7);
TH2F* hist32 = new TH2F("hist32","",100,600.0,1700.0,ybin*10,0,100);
  t1->Project("hist32",
   Form("(D.TH02-(%lf/sqrt(D.A02)+%lf))*0.05:D.A02",p22,q22),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
 
  hist32->GetXaxis()->SetRangeUser(1200,1600);
  hist32->GetYaxis()->SetRangeUser(30,40);
  hist32->GetXaxis()->SetTitle("ADC2 [ch]");
  hist32->GetYaxis()->SetTitle("TDC2 [ch/50 ps]");
  hist32->Draw("colz");
  c3->cd(8);



TH2F* hist33 = new TH2F("hist33","",100,600.0,1500.0,ybin*10,0,100);
  t1->Project("hist33",
   Form("(D.TH03-(%lf/sqrt(D.A03)+%lf))*0.05:D.A03",p23,q23),
	      cut_tdc00 && cut_tdc01 && cut_tdc02 && cut_tdc03 
	      && cut_qdc00 && cut_qdc01 && cut_qdc02 && cut_qdc03 
	      && cut_other);
 
  hist33->GetXaxis()->SetRangeUser(1100,1300);
  hist33->GetYaxis()->SetRangeUser(30,40);
  hist33->GetXaxis()->SetTitle("ADC3 [ch]");
  hist33->GetYaxis()->SetTitle("TDC3 [ch/50 ps]");
  hist33->Draw("colz");
  c0->cd(6);





  /////////////////////////////////////////////////////////////////////////////////////////


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
	      Form("((D.TH00+D.TH01)/2-(D.TH02+D.TH03)/2)*0.05-(%lf/sqrt(D.A00)+%lf)/2-(%lf/sqrt(D.A01)+%lf)/2+(%lf/sqrt(D.A02)+%lf)/2+(%lf/sqrt(D.A03)+%lf)/2",p20,q20,p21,q21,p22,q22,p23,q23),
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
  
 TF1* ffunc5 = new TF1("ffunc5","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-1,2.0);
 ffunc5->SetNpx(1000);
 ffunc5->SetParameters(1000,0,0.2);
 hist5->Fit("ffunc5","","",-1,1.5);
 hist5->GetXaxis()->SetRangeUser(-1,5);
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
hist20->Draw("colz");
  hist20->GetXaxis()->SetTitle("ADC0 [ch]");
  hist20->GetYaxis()->SetTitle("TDC0 [ch/50 ps]");
  hist20->GetXaxis()->SetRangeUser(900,1200);
  hist20->GetYaxis()->SetRangeUser(40,50);
  hist20->SetTitleSize(0.05,"x");
  hist20->SetTitleSize(0.05,"y");
  hist20->GetXaxis()->CenterTitle();
  hist20->GetYaxis()->CenterTitle();
hist20->Draw("colz");


  pf20->Draw("same");
  ffunc20->Draw("same");
 c2->cd(2);
 
hist30->Draw("colz");
  hist30->GetXaxis()->SetRangeUser(900,1200);
   hist30->GetYaxis()->SetRangeUser(40,50);
 hist30->GetXaxis()->SetTitle("ADC0 [ch]");
  
hist30->SetTitleSize(0.05,"x");
  hist30->SetTitleSize(0.05,"y");
hist30->GetYaxis()->SetTitle("TDC0 [ch/50 ps]");
  hist30->GetXaxis()->CenterTitle();
  hist30->GetYaxis()->CenterTitle(); 
hist30->Draw("colz");

 
}









