
void beta_ACnpe(){



  double beta,p,m;
  p=1200.; // particle momentum [MeV]
  m=934.8; // particle mass [MeV]
  beta=p/sqrt(pow(p,2)+pow(m,2)); // particle central beta
  cout<<"beta is "<<beta<<endl;
  int rnum;
  cout<<"run number : "<<endl;
  cin>>rnum;

  //======= Missing mass =============//
  double mm;// missing 

  /*
mm=sqrt(pow(,2)+pow(-L.tr.beta-R.tr.beta,2));
   */

  // ===================================== //                                               
  // ====== Open a ROOT file ============= //                                               
  // ===================================== //                                                
    TFile* f1 = new TFile(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d.root",rnum)); // position was displaed from center
  TTree* t1 = (TTree*)f1->Get("T");

  TCanvas* c0=new TCanvas("c0","c0");
  c0->Divide(3,1);
  double min_beta,max_beta,min_npe,max_npe;
  min_beta=-10.;
  max_beta=10.;
  min_npe=0.0;
  max_npe=1000.;
  int bin_beta,bin_npe;
  bin_npe=1000;
  bin_beta=20;
  TH2F*ha1=new TH2F("ha1","ha1",bin_npe,min_npe,max_npe,bin_beta,min_beta,max_beta);
  TH2F*ha2=new TH2F("ha2","ha2",bin_npe,min_npe,max_npe,bin_beta,min_beta,max_beta);

  TH1F* hbeta=new TH1F("hbeta","hbeta",bin_beta,min_beta,max_beta);
  t1->Project("ha1","R.tr.beta:R.a1.asum_c"
);
  c0->cd(1);
  ha1->Draw();
  c0->cd(2);
  t1->Project("ha2","R.tr.beta:R.a2.asum_c");
  ha2->Draw();
  c0->cd(3); 
 t1->Project("hbeta","R.tr.beta");
  hbeta->Draw();



}
