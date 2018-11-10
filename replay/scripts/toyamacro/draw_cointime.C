#include "Tree.cc"
#include "ParamMan.cc"
#include "Setting.cc"
#include "define.h"
void draw_cointime(){

  //TFile *ifp = new TFile("pdf/f1twc_check.root");
  TFile *ifp = new TFile("pdf/tmp.root");
  TTree *tree = (TTree*)ifp->Get("tree");
  Setting *set = new Setting();
  //ParamMan *param = new ParamMan("param/tmp.param");
  ParamMan *param = new ParamMan("param/f1_tuned_Lambda_twc.param");
  param->SetVal();
  int ENum = tree->GetEntries();

  TH2F *h2_ctime_lseg, *h2_ctime_rseg;
  TH1F *h_ctime_l[16], *h_ctime_r[16], *h_ctime;

  TF1 *ga_l[16], *ga_r[16];

  h2_ctime_lseg = new TH2F("h2_ctime_lseg","h2_ctime_lseg",1000,-10,10,16,-0.5,15.5);
  h2_ctime_rseg = new TH2F("h2_ctime_rseg","h2_ctime_rseg",1000,-10,10,16,-0.5,15.5);
  h_ctime = new TH1F("h_ctime","h_ctime",1000,-10,10);

  for(int i=0;i<16;i++){
    h_ctime_l[i] = new TH1F(Form("h_ctime_l%d",i),Form("h_ctime_l%d",i),1000,-10,10);
    h_ctime_r[i] = new TH1F(Form("h_ctime_r%d",i),Form("h_ctime_r%d",i),1000,-10,10);

  }

  double ctime;

  int s2lseg;
  double s2ltime,s2lat,s2lab;
  double s2ltt,s2ltb;
  double ltof;

  int s2rseg;
  double s2rtime,s2rat,s2rab;
  double s2rtt,s2rtb;
  double rtof;

  tree->SetBranchAddress("ctime"   ,&ctime  );
  tree->SetBranchAddress("s2lseg"  ,&s2lseg );
  tree->SetBranchAddress("s2ltime" ,&s2ltime);
  tree->SetBranchAddress("s2ltt"   ,&s2ltt  );
  tree->SetBranchAddress("s2ltb"   ,&s2ltb  );
  tree->SetBranchAddress("s2lat"   ,&s2lat  );
  tree->SetBranchAddress("s2lab"   ,&s2lab  );
  tree->SetBranchAddress("ltof"    ,&ltof   );
  tree->SetBranchAddress("s2rseg"  ,&s2rseg );
  tree->SetBranchAddress("s2rtime" ,&s2rtime);
  tree->SetBranchAddress("s2rtt"   ,&s2rtt  );
  tree->SetBranchAddress("s2rtb"   ,&s2rtb  );
  tree->SetBranchAddress("s2rat"   ,&s2rat  );
  tree->SetBranchAddress("s2rab"   ,&s2rab  );
  tree->SetBranchAddress("rtof"    ,&rtof   );

  for(int n=0;n<ENum;n++){
    if(n%10000==0)cout<<n <<" / "<<ENum<<endl;
    tree->GetEntry(n);
    h_ctime ->Fill(ctime);
    h2_ctime_lseg ->Fill(ctime,s2lseg);
    h2_ctime_rseg ->Fill(ctime,s2rseg);
    h_ctime_l[s2lseg]->Fill(ctime);
    h_ctime_r[s2rseg]->Fill(ctime);

  }

 int LR;
 for(int i=0;i<16;i++){
    //LR = 0;
    //ga_l[i] = new TF1(Form("ga_l%d",i+1),"gaus",-10,10);
    //set->SetTF1(ga_l[i],2,1,1);
    //double min=-50,max=50;
    //min = h_ctime_l[i]->GetXaxis()->GetBinCenter(h_ctime_l[i]->GetMaximumBin()) -3.;
    //max = h_ctime_l[i]->GetXaxis()->GetBinCenter(h_ctime_l[i]->GetMaximumBin()) +3.;
    //set->FitGaus(h_ctime_l[i],min,max,1.0,5);
    //h_ctime_l[i]->Fit(ga_l[i],"QR","",min,max);

    //param->SetTimeTune(CID_F1S2,i,LR,0,ga_l[i]->GetParameter(1) - 3.);
    //param->SetTimeTune(CID_F1S2,i,LR,1,ga_l[i]->GetParameter(1) - 3.);

    LR = 1;
    ga_r[i] = new TF1(Form("ga_r%d",i+1),"gaus",-10,10);
    set->SetTF1(ga_r[i],2,1,1);
    double min=-50,max=50;
    min = h_ctime_r[i]->GetXaxis()->GetBinCenter(h_ctime_r[i]->GetMaximumBin()) -3.;
    max = h_ctime_r[i]->GetXaxis()->GetBinCenter(h_ctime_r[i]->GetMaximumBin()) +3.;
    set->FitGaus(h_ctime_r[i],min,max,1.0,5);
    h_ctime_r[i]->Fit(ga_r[i],"QR","",min,max);

    param->SetTimeTune(CID_F1S2,i,LR,0,-1.*ga_r[i]->GetParameter(1) + 3.);
    param->SetTimeTune(CID_F1S2,i,LR,1,-1.*ga_r[i]->GetParameter(1) + 3.);
  }

  TCanvas *ca[3];
  for(int i=0;i<3;i++){
    ca[i] = new TCanvas(Form("ca%d",i+1), Form("c%d",i+1), 800, 800);
  }

  ca[0]->Clear(0);
  ca[0]->cd(1);
  h_ctime->Draw();

  ca[1]->Clear(0);
  ca[1]->cd(1);
  gPad->SetLogz(1);
  h2_ctime_lseg->Draw("colz");

  ca[2]->Clear(0);
  ca[2]->cd(1);
  gPad->SetLogz(1);
  h2_ctime_rseg->Draw("colz");

 //param->WriteToFile("param/tmp.param");
}
/*
  Tree *tr = new Tree();
  tr->add_tree("my_replay/tritium_111167.root");
  //tr->add_tree("my_replay/tritium_111168.root");
  //tr->add_tree("onl_replay_root/tritium_online_111145.root");
  tr->pack_tree();
  tr->readtreeS2L();
  tr->readtreeS2R();
  
  TH1F *h_coin = new TH1F("h_coin","h_coin",1000,0,300);

  int EN=tr->ENum;
  for(int n=0;n<EN;n++){
    tr->tree->GetEntry(n);

    for(int i=0;i<16;i++){
      for(int j=0;j<16;j++){
        if(tr->L_s2_lt[i]>0. &&tr->L_s2_rt[i]>0. &&tr->R_s2_lt[j]>0. &&tr->R_s2_rt[j]>0.){h_coin ->Fill(s2ns*tr->L_s2_time[i] - s2ns*tr->R_s2_time[j]);}//cout<< tr->L_s2_time[i] - tr->R_s2_time[i] <<endl;
      }
    }
  }

  h_coin->Draw();
*/
