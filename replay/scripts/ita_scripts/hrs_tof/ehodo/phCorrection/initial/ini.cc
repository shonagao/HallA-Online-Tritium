/*
  "ini.cc"
  
  Toshi Gogami , 15Oct2014
*/

#include<fstream>

void ini(){
  // ===== General conditions ========
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // ===== Open ROOT file ============
  TFile* f1 = new TFile("new_20171119.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  
  // ====== Cut conditions ============
  TCut cut1 = "adc1>500 && adc1<2000.0 && abs(time1-31.6)<3.0";
  TCut cut2 = "adc2>500 && adc2<2000.0 && abs(time2-21.3)<3.0";
  TCut cut3 = "adc3>500 && adc3<2000.0 && abs(time3-24.6)<3.0";
  TCut cut4 = "adc4>500 && adc4<2000.0 && abs(time4-22.8)<3.0";
  
  // ==== Create histograms =========
  TH2F* hist = new TH2F("hist","",
			100,500.0,1100.0,
			100,10.0,40.0);
  //hist->GetYaxis()->SetRangeUser(64.5,70.5);
  t1->Project("hist","time1:adc1",cut1&&cut2&&cut3&&cut4);
  hist->GetXaxis()->SetTitle("ADC [ch]");
  hist->GetYaxis()->SetTitle("TDC [ns]");
  
  TH2F* hist2 = (TH2F*)hist->Clone("hist2");
  //hist2->GetYaxis()->SetRangeUser(71.5,77.5);
  t1->Project("hist2","time2:adc2",cut1&&cut2&&cut3&&cut4);
  
  TH2F* hist3 = (TH2F*)hist->Clone("hist3");
  //hist3->GetYaxis()->SetRangeUser(68.0,74.0);
  t1->Project("hist3","time3:adc3",cut1&&cut2&&cut3&&cut4);
  
  TH2F* hist4 = (TH2F*)hist->Clone("hist4");
  //hist4->GetYaxis()->SetRangeUser(64.0,70.0);
  t1->Project("hist4","time4:adc4",cut1&&cut2&&cut3&&cut4);
  
  // ==== Profile =======
  TProfile* pf1 = new TProfile("pf1","",100,500.0,1400.0);
  t1->Project("pf1","time1:adc1",cut1&&cut2&&cut3&&cut4);
  pf1->SetLineColor(6);
  
  TProfile* pf2 = (TProfile*)pf1->Clone("pf2");
  t1->Project("pf2","time2:adc2",cut1&&cut2&&cut3&&cut4);
  
  TProfile* pf3 = (TProfile*)pf1->Clone("pf3");
  t1->Project("pf3","time3:adc3",cut1&&cut2&&cut3&&cut4);
  
  TProfile* pf4 = (TProfile*)pf1->Clone("pf4");
  t1->Project("pf4","time4:adc4",cut1&&cut2&&cut3&&cut4);
  
  // ==== Draw histograms ==========
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1); hist->Draw(); pf1->Draw("same");
  c1->cd(2); hist2->Draw(); pf2->Draw("same");
  c1->cd(3); hist3->Draw(); pf3->Draw("same");
  c1->cd(4); hist4->Draw(); pf4->Draw("same");
  
  
  // ===== Fitting ======
  c1->cd(1);
  //TF1* func1 = new TF1("func1","[0]+[1]*x+[2]*pow(x,2.0)",500.0,2500.0);
  //TF1* func1 = new TF1("func1","[0]+[1]*sqrt(x)+[2]*x ",500.0,2500.0);
  //TF1* func1 = new TF1("func1","[0]+[1]*sqrt(x)+[2]*sqrt(x*x) ",500.0,2500.0);
  TF1* func1 = new TF1("func1","[0]+[1]/sqrt(x)+[2]/x",500.0,2500.0);
  func1->SetParameters(50.0,1.0,1.0);
  pf1->Fit("func1","NB","",500.0,1200.0);
  func1->SetLineColor(9);
  func1->SetLineStyle(1);
  func1->Draw("same");
  c1->cd(2);
  TF1* func2 = (TF1*)func1->Clone("func2");
  func2->SetParameters(50.0,1.0,1.0);
  pf2->Fit("func2","NB","",500.0,1200.0);
  func2->Draw("same");
  c1->cd(3);
  TF1* func3 = (TF1*)func1->Clone("func3");
  func3->SetParameters(50.0,1.0,1.0);
  pf3->Fit("func3","NB","",500.0,1200.0);
  func3->Draw("same");
  c1->cd(4);
  TF1* func4 = (TF1*)func1->Clone("func4");
  func4->SetParameters(50.0,1.0,1.0);
  pf4->Fit("func4","NB","",500.0,1200.0);
  func4->Draw("same");
  
  const int n1 = 3; 
  double par1[n1];
  double par1er[n1];
  const int n2 = n1; 
  double par2[n2];
  double par2er[n2];
  const int n3 = n1; 
  double par3[n3];
  double par3er[n3];
  const int n4 = n1; 
  double par4[n4];
  double par4er[n4];
  for(int i=0 ; i<n1 ; i++){
    par1[i]   = func1->GetParameter(i);
    par1er[i] = func1->GetParError(i);
    par2[i]   = func2->GetParameter(i);
    par2er[i] = func2->GetParError(i);
    par3[i]   = func3->GetParameter(i);
    par3er[i] = func3->GetParError(i);
    par4[i]   = func4->GetParameter(i);
    par4er[i] = func4->GetParError(i);
  }
  
  // ===== Print ============
  //c1->Print("ini_toftest.eps","eps");
  //c1->Print("ini_toftest.png","png");
  
  // ===== Comments =======
  cout << endl;
  cout << " CH1: " << func1->GetExpFormula() << " " 
       << func1->GetNpar() << " "
       << par1[0] << " " << par1[1] << " " << par1[2] << endl;
  cout << " CH2: " << func2->GetExpFormula() << " " 
       << func2->GetNpar() << " "
       << par2[0] << " " << par2[1] << " " << par2[2] << endl;
  cout << " CH3: " << func3->GetExpFormula() << " " 
       << func3->GetNpar() << " "
       << par3[0] << " " << par3[1] << " " << par3[2] << endl;
  cout << " CH4: " << func4->GetExpFormula() << " " 
       << func4->GetNpar() << " "
       << par4[0] << " " << par4[1] << " " << par4[2] << endl;
  
  // ====== Output fitting results ==========
  ofstream* ofs = new ofstream("ini.dat");
  *ofs << " CH1: " << func1->GetExpFormula() << " " 
       << func1->GetNpar() << " "
       << par1[0] << " " << par1[1] << " " << par1[2] << endl;
  *ofs << " CH2: " << func2->GetExpFormula() << " " 
       << func2->GetNpar() << " "
       << par2[0] << " " << par2[1] << " " << par2[2] << endl;
  *ofs << " CH3: " << func3->GetExpFormula() << " " 
       << func3->GetNpar() << " "
       << par3[0] << " " << par3[1] << " " << par3[2] << endl;
  *ofs << " CH4: " << func4->GetExpFormula() << " " 
       << func4->GetNpar() << " "
       << par4[0] << " " << par4[1] << " " << par4[2] << endl;
  *ofs << endl;
  *ofs << " Initial parameters for pulse height correction." << endl;
  *ofs << " dragonvaio:~/ehodo/phCorrection/initial/ini.cc" << endl;
  *ofs << " Toshiyuki Gogami , 19Nov2017" << endl;
  ofs->close();
  
  cout << " --> Fitting results were output into: ini.dat" << endl;
  
}
