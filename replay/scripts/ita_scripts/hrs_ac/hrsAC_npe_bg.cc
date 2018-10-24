void hrsAC_npe(){

  gStyle->SetOptLogy(1);
  TFile* f1= new TFile("/adaqfs/home/a-onl/tritium_work/itabashi/ita_macro/rootfiles/tritium_new93403.root");
  TTree* t1 =(TTree*)f1->Get("tree");

  TH1F *hac1a[100]; //AC1 Ch24 Right &Left ADC
  TH1F *hac2a[100]; // AC2 ch26 Right & Left ADC
  TF1 *fac1a[100][2]; //AC1 Fitting func.
  TF1 *fac2a[100][2]; //AC2 Fitting func.

  double amin,amax;
  amin=-200.0;
  amax=1000.;
  int abin;
  abin=1000;
  double fac1min[24][2];
  double fac1max[24][2];
  double fac2min[26][2];
  double fac2max[26][2]; 

  //============== Seting Fitting Parameters ============== //
  //---- Initial parameters ------//
  //------ Setting Parameters ------------//

 //-------- AC1 ---------------//
  for(int i=0;i<24;i++){

    hac1a[i]=new TH1F(Form("hac1a[%d]",i),Form("hac1a[%d]",i),abin,amin,amax);
    t1->Project(Form("hac1a[%d]",i),Form("Ra1a_c[%d]",i)); //AC1 Right PMT

  //------Fit func. ----------//
    fac1a[i][0]=new TF1(Form("fac1a[%d]",i),gaus,From("fac1min[&d][0]",i),From("fac1max[%d][0]",i)); //pedestal
    fac1a[i][1]=new TF1(Form("fac1a[%d][1]",i),gaus,From("fac1min[&d][1]",i),From("fac1max[%d][1]",i));//1PE

    //   hac1a[i]->Fit(Form("fac1a[%d][0]",i));
    //  hac1a[i]->Fit(Form("fac1a[%d][1]",i));

}
  //------ AC2 ----------------//
  for(int i=0;i<26;i++){
    hac2a[i]=new TH1F(Form("hac2a[%d][0]",i),Form("hac2a[%d][0]",i),abin,amin,amax);
    t1->Project(Form("hac2a[%d]",i),Form("Ra2a_c[%d]",i)); //AC1
   
  //------ Fit func. ----------//    
    fac2a[i][0]=new TF1(Form("fac2a[%d][0]",i),gaus,Form("fac2min[%d][0]",i),Form("fmax[%d][0]",i)); //pedestal
    fac2a[i][1]=new TF1(Form("fac2a[%d][1]",i),gaus,Form("fac2min[%d][1]",i),Form("fac2max[%d][1]",i)); //1PE

    //  hac2a[i]->Fit(Form("fac2a[%d][0]",i));
    //  hac2a[i]->Fit(Form("fac2a[%d][1]",i));


}
 }


  //=========== Fitting ===================//



  TCanvas* c0 =new TCanvas("c0","c0");
  c0->Divide(2,2);
  c0->cd(1);
  hac1a[0][0]->Draw();
  c0->cd(2); 
  hac1a[0][1]->Draw();
  c0->cd(3);
  hac2a[0][0]->Draw();
  c0->cd(4); 
  hac2a[0][1]->Draw();




}
