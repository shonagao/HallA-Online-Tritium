
void sieveslit(int runnum=0){

  //-----------Calling rootfile & Settings-----------
  
  if(runnum==0){
    cout<<"Enter run number:"<<endl;
    cin>>runnum;
  }
  
  // TFile *fp = new TFile(Form("/adaqfs/home/a-onl/tritium_work/uehara/tohoku_analysis/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_online_%d.root",runnum));
  TFile *fp = new TFile(Form("/chafs1/work1/tritium/Rootfiles/tritium_online_%d.root",runnum));
  TTree *tree = (TTree*)fp->Get("T");
  
  int h_z_bin=400;
  int h_s_bin_th=400;
  int h_s_bin_ph=400;

  double h_z_start=-20;
  double h_z_end=20;
  double h_s_start_th=-0.08;
  double h_s_end_th=0.08;
  double h_s_start_ph=-0.05;
  double h_s_end_ph=0.05;
  double x_range=10;

  //-----------Making histgrams-----------

  TH1D *h_z = new TH1D("h_z","h_z",h_z_bin,h_z_start,h_z_end);
  SetTH1(h_z,"rpl.z","Z on beam line [cm]","Counts",1.,3002,2);

  // h_z->SetTitle("rpl.z");
  // h_z->SetLineColor(1);
  // h_z->SetLineWidth(1);
  // h_z->SetTitleSize(0.04,"");
  // h_z->SetTitleFont(42,"");
  // h_z->SetFillStyle(3002);
  // h_z->SetFillColor(2);

  // h_z->GetXaxis()->SetTitle("Z on beam line [cm]");
  // h_z->GetXaxis()->CenterTitle();
  // h_z->GetXaxis()->SetTitleFont(42);
  // h_z->GetXaxis()->SetTitleOffset(0.90);
  // h_z->GetXaxis()->SetTitleSize(0.06);
  // h_z->GetXaxis()->SetLabelFont(42);
  // h_z->GetXaxis()->SetLabelOffset(0.01);

  // h_z->GetYaxis()->SetTitle("Counts");
  // h_z->GetYaxis()->CenterTitle();
  // h_z->GetYaxis()->SetTitleFont(42);
  // h_z->GetYaxis()->SetTitleOffset(1.00);
  // h_z->GetYaxis()->SetTitleSize(0.06);
  // h_z->GetYaxis()->SetLabelFont(42);
  // h_z->GetYaxis()->SetLabelOffset(0.01);
  // ((TGaxis*)h_z->GetYaxis())->SetMaxDigits(4);

  tree->Project("h_z","rpl.z*100","");
  double x_max = h_z->GetXaxis()->GetBinCenter(h_z->GetMaximumBin());


  TH2D *h_s = new TH2D("h_s","h_s",h_s_bin_ph,h_s_start_ph,h_s_end_ph,h_s_bin_ph,h_s_start_th,h_s_end_th);
  SetTH2(h_s,"exL.th :exL.ph","Phi","Theta",0.);

  // h_s->SetTitle("sieveslit hole");
  // h_s->SetMinimum(0.8);
  // h_s->SetLineWidth(1);
  // h_s->SetTitleSize(0.05,"");
  // h_s->SetMarkerStyle(20);
  // h_s->SetMarkerSize(1.5);
  // h_s->SetMarkerColor(1);

  // h_s->GetXaxis()->SetTitle("Phi");
  // h_s->GetXaxis()->CenterTitle();
  // h_s->GetXaxis()->SetTitleFont(42);
  // h_s->GetXaxis()->SetTitleOffset(0.90);
  // h_s->GetXaxis()->SetTitleSize(0.06);
  // h_s->GetXaxis()->SetLabelFont(42);
  // h_s->GetXaxis()->SetLabelOffset(0.01);

  // h_s->GetYaxis()->SetTitle("Theta");
  // h_s->GetYaxis()->CenterTitle();
  // h_s->GetYaxis()->SetTitleFont(42);
  // h_s->GetYaxis()->SetTitleOffset(1.0);
  // h_s->GetYaxis()->SetTitleSize(0.06);
  // h_s->GetYaxis()->SetLabelFont(42);
  // h_s->GetYaxis()->SetLabelOffset(0.01);
  // ((TGaxis*)h_s->GetYaxis())->SetMaxDigits(4);

  tree->Project("h_s","exL.th : exL.ph","","");
  double h_s_max = h_s->GetBinContent(h_s->GetMaximumBin());

/*Form("rpl.z<%lf && rpl.z>%lf",x_max+x_range,x_max-x_range)*/
  //-----------Drawing-----------

  TCanvas *c = new TCanvas("c","c",1600,600);
  c->Divide(2,1,1E-5,1E-5);
  c->cd(1)->SetMargin(0.15,0.10,0.15,0.07);
  h_z->Draw();
  c->cd(2)->SetMargin(0.15,0.10,0.15,0.07);
  c->cd(2)->SetLogz();
  h_s->SetMinimum(1);
  // h_s->SetMaximum(h_s_max*0.5);
  // h_s->SetMarkerSize(0.01);
  h_s->Draw("colz");
}
