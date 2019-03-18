// 09/26/2018
// Since RHRS timing is too much in comparision with LHRS
// trying to see RHRRS RF-S2 time for the second time.
// Following exactely same steps from LHRS ie pad_7.c 
// Author Bishnu 
// Modify Itabashi 09/27/2018
// Pad# 7

void  hrs_tof_correct_posi()
{
  gStyle->SetOptStat(000000);
 TChain *T = new TChain("T");
 // HRS Right arm short range exp. run 
 for(Int_t i = 93495; i<93497;i++) /// these run goes 93295 to 93504
 {
   // T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

   T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));

 }

 //----------- Comment out --------------------------------//
 cout<<"============= < Right arm TOF posi analysis > ==================="<<endl;
 // cout<<"Channel number: "<<endl;
 // cin>>ch; 
 // cout<<"correction condition (x->0),(y->1),(theta->2),(phi->3):"<<endl;
 //cin>>tr;

 //------------------------------------------------------//
 double ch,chs;// convert ch->ns
 ch=0.0625;// [ns/ch]
 //chs=0.0625e-9;
 double xmin,xmax;
 int xbin;
 xmin=-1940.*ch;
 xmax=-1840.*ch;
 xbin=200;
 double fmin,fmax;
 fmin=-1905.*ch;
 fmax=-1885.*ch;
 int chmax =16;
 double p1[chmax],p2[chmax];
 double p1_err[chmax],p2_err[chmax];
 TH1F *htof[chmax]; 
 TF1 *ffit[chmax];

 for(int i=0;i<chmax;i++){

   htof[i]=new TH1F(Form("htof[%d]",i),Form("TOF (S2[%d]-S0)",i),xbin,xmin,xmax);
   T->Project(Form("htof[%d]",i),
	      Form("((R.s2.rt[%d]+R.s2.lt[%d])/2.0-(R.s0.rt+R.s0.lt)/2.0)*%lf",i,i,ch));
   ffit[i]=new TF1(Form("ffit[%d]",i),"landau",fmin,fmax); 
 
  //==== laudau function =====//
   /*
   ffit[i]=new TF1(Form("ffit[%d]",i),"landau",fmin,fmax); 
   htof[i]->Fit(Form("ffit[%d]",i),"","");
   p[i]=ffit[i]->GetParameter(1);
   */
  //==== gaussian function =====//
 
 ffit[i]=new TF1(Form("ffit[%d]",i),"gaus",fmin,fmax);  
  htof[i]->Fit(Form("ffit[%d]",i),"0","0");
   p1[i]=ffit[i]->GetParameter(1);
   p2[i]=ffit[i]->GetParameter(2);
   p1_err[i]=ffit[i]->GetParError(1);
   p2_err[i]=ffit[i]->GetParError(2);
 // htof[i]->Draw();
   htof[i]->SetTitle(Form("TOF times (S2[%d]-S0)",i));
   htof[i]->SetXTitle(Form("TOF(S2[%d]-S0) [ns]",i));
   htof[i]->SetYTitle("Counts/62.5 ps");
 }


 cout<<"======== TOF (S2-S0) infomation  ==========="<<endl;
 cout<<" TOF posi (S2[0]-S0): "<<p1[0]<<" [ns]"<<endl;
 cout<<"======== TOF position difference   ========="<<endl;
 for(int j=0;j<chmax;j++){
   cout<<"=======S2 channel ["<<j<<"] ============"<<endl;
   cout<<"TOF mean error :"<<p1_err[j]*1000<<"[ps]"<<endl;
   cout<<"TOF["<<j<<"]-TOF[0]"<<(p1[j]-p1[0])*1000<<" [ps]"<<endl;
   // cout<<"TOF sigma "<<p2[j]*1000<<" +- "<<p2_err[j]*1000<<" [ps]"<<endl; 
}
 double tof_max,tof_min;
 int kmin,kmax;
 tof_min=p1[0];
 tof_max=p1[0];
 for(int k=1;k<chmax;k++){
   if(p1[k]>tof_max){
     tof_max=p1[k];
     kmax=k;
   }
if(p1[k]<tof_min){
       tof_min=p1[k];
       kmin=k;
     }

 }

 TCanvas*c0 =new TCanvas("c0","c0");
 htof[kmin]->SetFillColor(kRed);
 htof[kmin]->Draw();
 htof[kmax]->SetFillColor(kBlue);
 htof[kmax]->Draw("same");

 cout<<"======== TOF most difference ========"<<endl;
 cout<<Form("TOF[%d]-TOF[%d] :",kmax,kmin)<<(tof_max-tof_min)*1000<<" [ps]"<<endl;

 }
