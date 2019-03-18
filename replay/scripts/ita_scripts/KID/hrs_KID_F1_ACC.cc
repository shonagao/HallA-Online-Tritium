

void hrs_KID_F1_ACC(){

  TChain*  T=new TChain("tree");
  tree->Add("/home/itabashi/jlab_nnL/ita_scripts/rootfiles/coin_H2_1.root");

  double ctime[100],sum_a1,sum_a2;
 T->SetBranchStatus("*",0);  
 T->SetBranchStatus("ctime",1);    
 T->SetBranchAddress("ctime",ctime);
 T->SetBranchStatus("R.a1.asum_c",1);    
 T->SetBranchAddress("R.a1.asum_c",&sum_a1);
 T->SetBranchStatus("R.a2.asum_c",1);    
 T->SetBranchAddress("R.a2.asum_c",&sum_a2);
 TH1F* hcoin_acc=new TH1F("hcoin_acc","Coin ACC region Hist",1000,-20,20.);


 int evnt=tree->GetEntries();
 cout<<"Event num :"<<evnt<<endl;
 for(int k=0;k<evnt;k++){
   tree->GetEntry(k);
   if(2.0<sum_a2 && sum_a2<7.0)hcoin_acc->Fill(ctime[0]);
}


 int max;
 max=3;
 TF1* facc[1000];
 for(int i=0;i<max;i++){
 facc[i]=new TF1(Form("facc[%d]",i),"[0]+gaus(1)",-16+2*i,-14+2*i);
 facc[i]->SetParameter(0,60.);
 facc[i]->SetParLimits(1,10,100);
 facc[i]->SetParameter(2,-16.+2*i);
 facc[i]->SetParLimits(3,0.1,0.6);
 hcoin_acc->Fit(Form("facc[%d]",i),"","",-17+2*i,-15+2*i);
 }


 TCanvas* c0=new TCanvas("c0","Coin time ACC region Hist");
 c0->cd();
 hcoin_acc->Draw();
 for(int i=0;i<max;i++){
 facc[i]->SetLineColor(i+1);
 facc[i]->Draw("same");
 }

}
