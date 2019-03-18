double mean_1PE(int AC,int SEG);

void hrsAC_npe_ana(){

    bool print_flag=true;
    if(print_flag)gROOT->SetBatch(0);
    gStyle->SetOptLogy(1);



  TChain*  T=new TChain("T");
  // T->Add("/data/7b/E12-17-003/replay/Rootfiles/nnL/tritium_111621.root"); //cosmic run
  T->Add("/data/7b/E12-17-003/replay/Rootfiles/nnL/tritium_111623.root"); //Tritium run
  double Ra1a_p[24],Ra2a_p[26],Ra1a_c[24],Ra2a_c[26],Ra1asum_p,Ra2asum_p;
 T->SetBranchStatus("*",0);  
 T->SetBranchStatus("R.a1.a_p",1);
 T->SetBranchAddress("R.a1.a_p",Ra1a_p);
 T->SetBranchStatus("R.a2.a_p",1);
 T->SetBranchAddress("R.a2.a_p",Ra2a_p);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c);
 T->SetBranchStatus("R.a1.asum_p",1);
 T->SetBranchAddress("R.a1.asum_p",Ra1asum_p);
 T->SetBranchStatus("R.a2.asum_p",1);
 T->SetBranchAddress("R.a2.asum_p",Ra2asum_p);
 TH1F* ha1a_p[24];
 TH1F* ha2a_p[26];
 TH1F* ha1a_c[24];
 TH1F* ha2a_c[26];
 TH1F* ha1a_npe[24];
 TH1F* ha2a_npe[26];
 
 TH1F* ha1_sum=new TH1F("ha1_sum","AC1 NPE SUM",500,0.0,100.0);
 TH1F* ha2_sum=new TH1F("ha2_sum","AC2 NPE SUM",500,0.0,100.0);

 double min_ac1_p,max_ac1_p,bin_ac1_p;
 min_ac1_p=-200.0;
 max_ac1_p=3000.;
 bin_ac1_p=max_ac1_p-min_ac1_p;
 bin_ac1_p=(int)bin_ac1_p;

 cout<<"bin_ac1_p"<<bin_ac1_p<<endl;
 double min_ac2_p,max_ac2_p,bin_ac2_p;
 min_ac2_p=-200.0;
 max_ac2_p=3000.;
 bin_ac2_p=max_ac2_p-min_ac2_p;
 bin_ac2_p=(int)bin_ac2_p;
 cout<<"bin_ac2_p"<<bin_ac2_p<<endl;
 double min_ac1_npe,max_ac1_npe,bin_ac1_npe;
 min_ac1_npe=-5.0;
 max_ac1_npe=50.;
 bin_ac1_npe=max_ac1_npe-min_ac1_npe;
 bin_ac1_npe=(int)bin_ac1_npe;
 double min_ac2_npe,max_ac2_npe,bin_ac2_npe;
 min_ac2_npe=-5.0;
 max_ac2_npe=100.;
 bin_ac2_npe=max_ac2_npe-min_ac2_npe;
 bin_ac2_npe=(int)bin_ac2_npe;

 double ac1_sum,ac2_sum;
 ac1_sum=0.0;
 ac2_sum=0.0;
 int evnt= T->GetEntries();
 cout<<"Get Entries: "<<evnt<<endl;
 int c=1;

 for(int i=0;i<24;i++){
     ha1a_p[i]=new TH1F(Form("ha1a_p[%d]",i),Form("AC1 ADC Hist Seg%d",i),bin_ac1_p,min_ac1_p,max_ac1_p);
     // ha1a_c[i]=new TH1F(Form("ha1a_c[%d]",i),Form("AC1 NPE Hist Seg%d",i),bin_ac1_c,min_ac1_c,max_ac1_c);
     ha1a_npe[i]=new TH1F(Form("ha1a_npe[%d]",i),Form("AC1 NPE Hist Seg%d",i),bin_ac1_npe,min_ac1_npe,max_ac1_npe);
}

   for(int j=0;j<26;j++){ 
     ha2a_p[j]=new TH1F(Form("ha2a_p[%d]",j),Form("AC2 ADC Hist Seg%d",j),bin_ac2_p,min_ac2_p,max_ac2_p);
     ha2a_npe[j]=new TH1F(Form("ha2a_npe[%d]",j),Form("AC2 NPE Hist Seg%d",j),bin_ac2_npe,min_ac2_npe,max_ac2_npe);
     // ha2a_c[j]=new TH1F(Form("ha2a_c[%d]",j),Form("AC2 NPE Hist Seg%d",j),bin_ac2_c,min_ac2_c,max_ac2_c);
   }


  for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   if(k==c*10000){
     cout<<"Event : "<<k<<"/ "<<evnt<<endl;
     c=c+1;
     }
   for(int i=0;i<24;i++){
     ha1a_p[i]->Fill(Ra1a_p[i]);
     ha1a_npe[i]->Fill(Ra1a_p[i]/mean_1PE(1,i));
     ac1_sum=ac1_sum+Ra1a_p[i]/mean_1PE(1,i);
     }
    for(int j=0;j<26;j++){
     ha2a_p[j]->Fill(Ra2a_p[j]);
     ha2a_npe[j]->Fill(Ra2a_p[j]/mean_1PE(2,j));
     ac2_sum=ac2_sum+Ra2a_p[j]/mean_1PE(2,j);
     }

    ha1_sum->Fill(ac1_sum);
    ha2_sum->Fill(ac2_sum);

    ac1_sum=0.0;
    ac2_sum=0.0;

}


 cout<<"Fill Dene !!"<<endl;

 TF1 *fa1_c[24];
 TF1 *fa2_c[26];
 TF1 *fa1_p[24];
 TF1 *fa2_p[26];
 TGraph* ga1=new TGraph(24);
 TGraph* ga2=new TGraph(26);
 double mean_a1_c[24],mean_a2_c[26],mean_a1_p[24],mean_a2_p[26]; 
 for(int i=0;i<24;i++){
 fa1_p[i]=new TF1(Form("fa1_p[%d]",i),"gaus",100.,600.);
 ha1a_p[i]->Fit(Form("fa1_p[%d]",i),"Rq","",100.,600.);
 mean_a1_p[i]=fa1_p[i]->GetParameter(2);
 ga1->SetPoint(i,i,mean_a1_p[i]);
 // cout<<"AC1 seg"<<i<<" 1PE Mean:"<<mean_a1_p[i]<<endl;
 }
 for(int j=0;j<26;j++){
 fa2_p[j]=new TF1(Form("fa2_p[%d]",j),"gaus",100.,600.);
 ha2a_p[j]->Fit(Form("fa2_p[%d]",j),"Rq","",100.,600.);
 mean_a2_p[j]=fa2_p[j]->GetParameter(2);
 ga2->SetPoint(j,j,mean_a2_p[j]);
 //cout<<"AC2 seg"<<j<<" 1PE Mean:"<<mean_a2_p[j]<<endl;
 }








 TCanvas* c0=new TCanvas("c0","ha1a_p hist");
 TCanvas* c1=new TCanvas("c1","ha1a_p hist2");
 c0->Divide(4,3);
 c1->Divide(4,3);

 for(int i=0;i<12;i++){
   c0->cd(i+1);
   ha1a_p[i]->Draw();
   fa1_p[i]->SetLineColor(2);
   fa1_p[i]->Draw("same");
   c1->cd(i+1);
   ha1a_p[i+12]->Draw();
   fa1_p[i+12]->SetLineColor(2);
   fa1_p[i+12]->Draw("same");
}
 

 TCanvas* c2=new TCanvas("c2","ha2a_p hist ");
 TCanvas* c3=new TCanvas("c3","ha2a_p hist2 ");
 c2->Divide(4,4);
 c3->Divide(4,4);
 for(int i=0;i<13;i++){
   c2->cd(i+1);
   ha2a_p[i]->Draw();
   fa2_p[i]->SetLineColor(2);
   fa2_p[i]->Draw("same");
  c3->cd(i+1);
   ha2a_p[i+13]->Draw();
   fa2_p[i+13]->SetLineColor(2);
   fa2_p[i+13]->Draw("same");
}

 TCanvas* c4 =new TCanvas("c4","AC 1PE mean");
 c4->Divide(1,2);
 c4->cd(1);
 c4->cd(1)->SetLogy(0);
 ga1->SetMarkerStyle(21);
 ga1->SetMarkerColor(kRed);
 ga1->SetMarkerSize(0.45);
 ga1->Draw("AP");
 c4->cd(2);
 c4->cd(2)->SetLogy(0);
 ga2->SetMarkerStyle(21);
 ga2->SetMarkerColor(kRed);
 ga2->SetMarkerSize(0.45);
 ga2->Draw("AP");

 TCanvas* c5=new TCanvas("c5","AC NPE SUM Hist");
 c5->Divide(1,2);
 c5->cd(1);
 ha1_sum->Draw();
 c5->cd(2);
 ha2_sum->Draw();

TCanvas* c6= new TCanvas("c6","AC1 NPE Hist");
 c6->Divide(4,3);
TCanvas* c7= new TCanvas("c7","AC1 NPE Hist2");
 c7->Divide(4,3);
 for(int i=0;i<12;i++){
   c6->cd(i+1);
   ha1a_npe[i]->Draw();
   c7->cd(i+1);
   ha1a_npe[i+12]->Draw();
}

TCanvas* c8= new TCanvas("c8","AC1 NPE Hist");
 c8->Divide(4,4);
TCanvas* c9= new TCanvas("c9","AC1 NPE Hist2");
 c9->Divide(4,4);
 for(int i=0;i<13;i++){
   c8->cd(i+1);
   ha2a_npe[i]->Draw();
   c9->cd(i+1);
   ha2a_npe[i+13]->Draw();
}

 cout<<"Drawn !!"<<endl;


 if(print_flag){
 string name="AC_npe_ana.pdf";
 c0->Print(Form("%s[",name.c_str()),"pdf");
 c0->Print(Form("%s",name.c_str()),"pdf");
 c1->Print(Form("%s",name.c_str()),"pdf");
 c2->Print(Form("%s",name.c_str()),"pdf");
 c3->Print(Form("%s",name.c_str()),"pdf");
 c4->Print(Form("%s",name.c_str()),"pdf");
 c5->Print(Form("%s",name.c_str()),"pdf");
 c6->Print(Form("%s",name.c_str()),"pdf");
 c7->Print(Form("%s",name.c_str()),"pdf");
 c8->Print(Form("%s",name.c_str()),"pdf");
 c9->Print(Form("%s",name.c_str()),"pdf");
 c9->Print(Form("%s]",name.c_str()),"pdf");
}

 cout<<"Finish making pdf "<<endl;

 //======== Comment Out =============//
 cout<<"AC1 1PE Mean [ch]:{ ";
 for(int i=0;i<24;i++){
   cout<<mean_a1_p[i]<<",";
}
 cout<<"}"<<endl;

 cout<<"AC2 1PE Mean [ch]:{ ";
 for(int i=0;i<26;i++){
   cout<<mean_a2_p[i]<<",";
}
 cout<<"}"<<endl;


 gSystem->Exit(1);



 /*
for(int i=0;i<24;i++){
 ha1a_c[i]->Fit(Form("fa1_c[%d]",i),"gaus",0.5.,1.5);
 mean_a1_c[i]=fa1_c[i]->GetParameter(3);
 cout<<"AC1 seg"<<i<<" 1PE Mean:"<<mean_a1_c[i]<<endl;
 }
 for(int j=0;j<26;j++){
 ha2a_c[j]->Fit(Form("fa2_c[%d]",i),"gaus",0.5.,1.5);
 mean_a2_c[j]=fa2_c[i]->GetParameter(3);
 cout<<"AC2 seg"<<j<<" 1PE Mean:"<<mean_a2_c[j]<<endl;
 }
 */




}



//===== Parameters =============//
double mean_1PE(int AC,int SEG){

  double adc_mean[26];
  if(AC==1){
double AC1_mean[24]={166.715,129.851,131.465,116.276,129.286,155.609,116.604,125.268,111.062,124.192,128.875,140.905,92.367,135.241,96.183,96.3441,110.618,115.393,125.49,143.868,136.808,122.388,131.726,136.755}
 adc_mean[SEG]=AC1_mean[SEG];
  }else if(AC==2){
    double AC2_mean[26]={205.764,174.776,181.969,184.082,170.038,206.409,208.37,165.576,221.58,172.785,203.037,187.834,156.645,170.033,150.977,152.156,184.658,202.388,190.307,240.03,202.344,175.646,201.902,171.968,129.016,225.013}

    adc_mean[SEG]=AC2_mean[SEG];
  }

  return adc_mean[SEG];

}
