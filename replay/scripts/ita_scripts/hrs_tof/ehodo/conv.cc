/*
  conv.cc
  
  Toshiyuki Gogami, November 18, 2017
*/

const double ch2ns = 0.05; // ch --> ns 

void conv(){
  // =======
  TFile* f1 = new TFile("hkstest_01330.root");
  TTree* t1 = (TTree*)f1->Get("T");

  const int nn = 50;
  double adc1, adc2, adc3, adc4, adc5, adc6;
  double tdc1[nn], tdc2[nn], tdc3[nn], tdc4[nn], tdc5[nn], tdc6[nn];
  double tof;
  //double tdc1, tdc2, tdc3, tdc4;
  int n1, n2, n3, n4, n5, n6;
  t1->SetBranchAddress("D.A00",&adc1);
  t1->SetBranchAddress("D.A01",&adc2);
  t1->SetBranchAddress("D.A02",&adc3);
  t1->SetBranchAddress("D.A03",&adc4);
  t1->SetBranchAddress("D.A04",&adc5);
  t1->SetBranchAddress("D.A05",&adc6);
  t1->SetBranchAddress("Ndata.D.TH00",&n1);
  t1->SetBranchAddress("Ndata.D.TH01",&n2);
  t1->SetBranchAddress("Ndata.D.TH02",&n3);
  t1->SetBranchAddress("Ndata.D.TH03",&n4);
  t1->SetBranchAddress("Ndata.D.TH04",&n5);
  t1->SetBranchAddress("Ndata.D.TH05",&n6);
  t1->SetBranchAddress("D.TH00",&tdc1);
  t1->SetBranchAddress("D.TH01",&tdc2);
  t1->SetBranchAddress("D.TH02",&tdc3);
  t1->SetBranchAddress("D.TH03",&tdc4);
  t1->SetBranchAddress("D.TH04",&tdc5);
  t1->SetBranchAddress("D.TH05",&tdc6);

  int evid = 0;
  TFile* fnew = new TFile("new.root","recreate");
  TTree* tnew = new TTree("tree","EHODO test");
  tnew->Branch("evid",&evid,"evid/I");
  tnew->Branch("adc1",&adc1,"adc1/D");
  tnew->Branch("adc2",&adc2,"adc2/D");
  tnew->Branch("adc3",&adc3,"adc3/D");
  tnew->Branch("adc4",&adc4,"adc4/D");
  tnew->Branch("adc5",&adc5,"adc5/D");
  tnew->Branch("adc6",&adc6,"adc6/D");
  //tnew->Branch("tdc1",&tdc1[0],"tdc1/D");
  //tnew->Branch("tdc2",&tdc2[0],"tdc2/D");
  //tnew->Branch("tdc3",&tdc3[0],"tdc3/D");
  //tnew->Branch("tdc4",&tdc4[0],"tdc4/D");
  //tnew->Branch("tdc5",&tdc5[0],"tdc5/D");
  //tnew->Branch("tdc6",&tdc6[0],"tdc6/D");
  tnew->Branch("time1",&tdc1[0],"time1/D");
  tnew->Branch("time2",&tdc2[0],"time2/D");
  tnew->Branch("time3",&tdc3[0],"time3/D");
  tnew->Branch("time4",&tdc4[0],"time4/D");
  tnew->Branch("time5",&tdc5[0],"time5/D");
  tnew->Branch("time6",&tdc6[0],"time6/D");
  tnew->Branch("tof",&tof,"tof/D");
  
  double ent = t1->GetEntries();
  for(int i=0 ; i<ent ; i++){
    evid = i;
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n4 = 0;
    n5 = 0;
    n6 = 0;
    adc1 = -2222.0;
    adc2 = -2222.0;
    adc3 = -2222.0;
    adc4 = -2222.0;
    adc5 = -2222.0;
    adc6 = -2222.0;
    for(int j=0; j<50; j++){
      tdc1[j] = -2222.0;
      tdc2[j] = -2222.0;
      tdc3[j] = -2222.0;
      tdc4[j] = -2222.0;
      tdc5[j] = -2222.0;
      tdc6[j] = -2222.0;
    }
    t1->GetEntry(i);
    for(int j=0 ; j<50 ; j++){
      tdc1[j] = tdc1[j] * ch2ns;
      tdc2[j] = tdc2[j] * ch2ns;
      tdc3[j] = tdc3[j] * ch2ns;
      tdc4[j] = tdc4[j] * ch2ns;
      tdc5[j] = tdc5[j] * ch2ns;
      tdc6[j] = tdc6[j] * ch2ns;
    }
    tof = (tdc1[0]+tdc2[0])/2.0 - (tdc3[0]+tdc4[0])/2.0;
    //cout << adc1 << endl;
    //cout << n1 << endl;
    //cout << tdc1[0] << endl;
    if(n1==1 && n2==1 && n3==1 && n4==1){
      tnew->Fill();
    }
  }
  tnew->Write();
  fnew->Close();
  
}
