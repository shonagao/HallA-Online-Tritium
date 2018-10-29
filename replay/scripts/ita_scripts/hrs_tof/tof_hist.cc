// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void tof_hist(){
  TChain *tree = new TChain("T");
// tree->Add("rootfiles/tritium_online_94000.root");
//  tree->Add("rootfiles/tritium_online_93926.root");
  int run;
  cout<<"run number : ";
  cin>>run;

  tree->Add(Form("/u/home/itabashi/nnL/Tohoku_github/HallA-Online-Tritium/replay/t2root/Rootfiles/tritium_%d.root",run));

  TCanvas *c1 = new TCanvas("c1","c1",800,600); c1->Divide(2,8);
  TCanvas *cbeta = new TCanvas("cbeta","cbeta",800,600); cbeta->Divide(2,8);

  TH1D *h1[16];
  TH1D *hbeta[16];
  string F1 = "RTDC.F1FirstHit";
  double fact = 0.05623;//[ns]
  string S0t = Form("%s[%d]",F1.c_str(),43);
  string S0b = Form("%s[%d]",F1.c_str(),44);
  string S2l[16];
  string S2r[16];
  string L1  = Form("%s[%d]",F1.c_str(),9);
  string L2  = Form("%s[%d]",F1.c_str(),46);
  string pathl = "R.tr.pathl";
  string cut[16];
  

  for(int i=0;i<16;i++){

  S2l[i]= Form("%s[%d]",F1.c_str(),16+i);
  S2r[i]= Form("%s[%d]",F1.c_str(),48+i);
  cut[i] = Form("%s>0 && %s>0 && %s>0 && %s>0 && %s>0 && %s>0"
            ,S0t.c_str(),S0b.c_str(),S2l[i].c_str(),S2r[i].c_str(),L1.c_str(),L2.c_str());


    h1[i] = new TH1D(Form("h1[%d]",i),Form("h1[%d]",i),1000,-300.,0);
    hbeta[i] = new TH1D(Form("hbeta[%d]",i),Form("hbeta[%d]",i),200,-1.,1);

  tree->Project(Form("h1[%d]",i),Form("(((%s-%s)+(%s-%s))/2. - ((%s-%s)+(%s-%s))/2.)*%lf"
                    ,S2l[i].c_str(),L1.c_str(),S2r[i].c_str(),L2.c_str()
                    ,S0t.c_str(),L2.c_str(),S0b.c_str(),L2.c_str(),fact)
                    ,Form("%s",cut[i].c_str()));

  tree->Project(Form("hbeta[%d]",i),Form("%s/3.0e8/((((%s-%s)+(%s-%s))/2. - ((%s-%s)+(%s-%s))/2.0)*%lf*1.0e-9)"
		    ,pathl.c_str(),S2l[i].c_str(),L1.c_str(),S2r[i].c_str(),L2.c_str()
                    ,S0t.c_str(),L2.c_str(),S0b.c_str(),L2.c_str(),fact)
	            ,Form("%s",cut[i].c_str()));

  c1->cd(i+1);
  gPad->SetLogy();
  h1[i]->SetXtitle(Form("S2ch%d-S0 Hist [ns]"),i);
  h1[i]->Draw();
  cbeta->cd(i+1);
  hbeta[i]->SetXtitle(Form("S2ch%d -S0 pathl beta Hist"),i);
  hbeta[i]->Draw();

  }

  TString name;
  name=Form("./pdf/tof_hist_%d.pdf",run);
  c1->Print(name +"[");
  c1->Print(name);
  cbeta->Print(name);
  cbeta->Print(name+"]");
    }

