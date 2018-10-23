/*
  adc_tof.cc
  Toshi Gogami , 13Oct2014
*/

void adc_tof(){
  // ====== General conditions =======
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1011);
  
  // ====== Open ROOT file ============
  TFile* f1 = new TFile("phc_test.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  
  // ====== Cut conditions ============
  TCut cut1 = "adc1>500 && adc1<2000.0 && abs(time1-31.6)<5.0";
  TCut cut2 = "adc2>500 && adc2<2000.0 && abs(time2-21.3)<5.0";
  TCut cut3 = "adc3>590 && adc3<2000.0 && abs(time3-24.6)<5.0";
  TCut cut4 = "adc4>600 && adc4<2000.0 && abs(time4-22.8)<5.0";
  
  // ====== Create histograms =========
  double amin = 0.0 , amax = 3000.0;
  //double tmin = 50.0 , tmax = 100.0; // time
  double tmin = -25.0 , tmax = 25.0; // ctime
  int abin = 375;
  int tbin = 100;
  // ~~~ ADC1 vs. TDC1~~~~
  TH2F* hist1 = new TH2F("hist1","",
			 abin,amin,amax,
			 tbin,tmin,tmax);
  //t1->Project("hist1","time1:adc1",cut1&&cut2&&cut3&&cut4);
  t1->Project("hist1","ctime1:adc1",cut1&&cut2&&cut3&&cut4);
  //t1->Project("hist1","time1:adc1",cut1&&cut2&&cut3);
  hist1->GetXaxis()->SetRangeUser(550.0,1000.);
  hist1->GetYaxis()->SetRangeUser(-5.0,5.0);
  hist1->GetXaxis()->SetTitle("ADC [ch]");
  hist1->GetYaxis()->SetTitle("Time [ns]");
  // ~~~ ADC2 vs. TDC2~~~~
  TH2F* hist2 = (TH2F*)hist1->Clone("hist2");
  //hist2->GetYaxis()->SetRangeUser(71.0,78.);
  //t1->Project("hist2","time2:adc2",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist2","ctime2:adc2",cut1&&cut2&&cut3&&cut4); 
  // ~~~ ADC3 vs. TDC3~~~~
  TH2F* hist3 = (TH2F*)hist1->Clone("hist3");
  //hist3->GetYaxis()->SetRangeUser(65.0,74.);
  //t1->Project("hist3","time3:adc3",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist3","ctime3:adc3",cut1&&cut2&&cut3&&cut4);  
  // ~~~ ADC4 vs. TDC4~~~~
  TH2F* hist4 = (TH2F*)hist1->Clone("hist4");
  //hist4->GetYaxis()->SetRangeUser(62.0,70.);
  //t1->Project("hist4","time4:adc4",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist4","ctime4:adc4",cut1&&cut2&&cut3&&cut4);  
  // ~~~ Time Of Flight ~~~~
  TH1F* h_tof = new TH1F("h_tof","TOF",100,-2.5,2.5);
  //h_tof->GetYaxis()->SetRangeUser(62.0,70.);
  //t1->Project("h_tof","tof-2.5",cut1&&cut2&&cut3&&cut4);
  t1->Project("h_tof","ctof",cut1&&cut2&&cut3&&cut4);
  h_tof->GetXaxis()->SetTitle("[ns]");
  //h_tof->GetYaxis()->SetTitle("[Counts / 100 ps]");
  h_tof->GetYaxis()->SetTitle("[Counts / 0.05 ns]");
  h_tof->SetFillColor(9);
  
  
  // ~~~ ADC1 vs. TOF~~~~
  //double tofmin = -5.0 , tofmax = 5.0;
  double tofmin = -2.5 , tofmax = 2.5;
  int tofbin = 100;
  TH2F* hist5 = new TH2F("hist5","",
			 abin,amin,amax,
			 tofbin,tofmin,tofmax);
  //t1->Project("hist5","tof-2.5:adc1",cut1&&cut2&&cut3&&cut4);
  t1->Project("hist5","ctof:adc1",cut1&&cut2&&cut3&&cut4);
  hist5->GetYaxis()->SetRangeUser(-1.5,1.5);
  hist5->GetXaxis()->SetRangeUser(550.0,1000.);
  hist5->GetXaxis()->SetTitle("ADC [ch]");
  hist5->GetYaxis()->SetTitle("TOF [ns]");
  // ~~~ ADC2 vs. TOF~~~~
  TH2F* hist6 = (TH2F*)hist5->Clone("hist6");
  //t1->Project("hist6","tof-2.5:adc2",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist6","ctof:adc2",cut1&&cut2&&cut3&&cut4);  
  // ~~~ ADC3 vs. TOF~~~~
  TH2F* hist7 = (TH2F*)hist5->Clone("hist7");
  //t1->Project("hist7","tof-2.5:adc3",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist7","ctof:adc3",cut1&&cut2&&cut3&&cut4);  
  // ~~~ ADC4 vs. TOF~~~~
  TH2F* hist8 = (TH2F*)hist5->Clone("hist8");
  //t1->Project("hist8","tof-2.5:adc4",cut1&&cut2&&cut3&&cut4);  
  t1->Project("hist8","ctof:adc4",cut1&&cut2&&cut3&&cut4);  
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);
  //t1->Draw("time1:adc1");
  c1->cd(1);hist1->Draw("colz");
  c1->cd(2);hist2->Draw("colz");
  c1->cd(3);hist3->Draw("colz");
  c1->cd(4);hist4->Draw("colz");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  h_tof->Draw();

  // ===== Fitting =======
  //TF1* func1 = new TF1("func1","gaus",-3.0,3.0);
  //TF1* func1 = new TF1("func1","[0]*exp(x)",-3.0,3.0);
  TF1* func1 = new TF1("func1","[0]*1.0/sqrt(2.0*3.14159)/[2]*exp(-0.5*pow((x-[1])/[2],2.0))",-2.0,2.0);
  func1->SetParameters(500.0,0.0,0.13);
  //func1->SetParLimits(2,0.12,0.16);
  //func1->FixParameter(1,0.);
  func1->SetNpx(1000);
  //h_tof->Fit("func1","","",-1.5.1.5);
  //h_tof->Fit("func1","","",-0.8.1.8);
  //h_tof->Fit("func1","","",-1.0.1.0);
  //h_tof->Fit("func1","BN","",-0.15.0.3); 
  //h_tof->Fit("func1","BN","",-3.0.3.0); 
  h_tof->Fit("func1","BN","",-0.3,0.3); 
  func1->SetLineColor(2);
  func1->SetLineWidth(2);
  //func1->SetLineStyle(7);
  func1->Draw("same");

  TCanvas* c3 = new TCanvas("c3","c3");
  c3->Divide(2,2);
  c3->cd(1);hist5->Draw("colz");
  c3->cd(2);hist6->Draw("colz");
  c3->cd(3);hist7->Draw("colz");
  c3->cd(4);hist8->Draw("colz");
  
  // ===== Print ======
  //c1->Print("pict/toftest_adcVStdc_1-1-1-1.eps","eps");
  //c1->Print("pict/toftest_adcVStdc_1-1-1-1.png","png");
  //c2->Print("pict/toftest_tof_1-1-1-1.eps","eps");
  //c2->Print("pict/toftest_tof_1-1-1-1.png","png");
  //c3->Print("pict/toftest_adcVStof_1-1-1-1.eps","eps");
  //c3->Print("pict/toftest_adcVStof_1-1-1-1.png","png");
  
  //c1->Print("pict/toftest_adcVStdc_0-1-1-1.eps","eps");
  //c1->Print("pict/toftest_adcVStdc_0-1-1-1.png","png");
  //c2->Print("pict/toftest_tof_0-1-1-1.eps","eps");
  //c2->Print("pict/toftest_tof_0-1-1-1.png","png");
  //c3->Print("pict/toftest_adcVStof_0-1-1-1.eps","eps");
  //c3->Print("pict/toftest_adcVStof_0-1-1-1.png","png");
  
  //c1->Print("pict/toftest_adcVStdc_0-1-1-4.eps","eps");
  //c1->Print("pict/toftest_adcVStdc_0-1-1-4.png","png");
  //c2->Print("pict/toftest_tof_0-1-1-4.eps","eps");
  //c2->Print("pict/toftest_tof_0-1-1-4.png","png");
  //c3->Print("pict/toftest_adcVStof_0-1-1-4.eps","eps");
  //c3->Print("pict/toftest_adcVStof_0-1-1-4.png","png");
  
  //c1->Print("pict_tofadc/toftest_adcVStdc_1-1-1-1-0.eps","eps");
  //c1->Print("pict_tofadc/toftest_adcVStdc_1-1-1-1-0.png","png");
  //c2->Print("pict_tofadc/toftest_tof_1-1-1-1-0.eps","eps");
  //c2->Print("pict_tofadc/toftest_tof_1-1-1-1-0.png","png");
  //c3->Print("pict_tofadc/toftest_adcVStof_1-1-1-1-0.eps","eps");
  //c3->Print("pict_tofadc/toftest_adcVStof_1-1-1-1-0.png","png");
}
