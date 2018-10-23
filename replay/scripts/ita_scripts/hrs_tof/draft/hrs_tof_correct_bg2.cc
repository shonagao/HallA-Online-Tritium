// 09/26/2018
// Since RHRS timing is too much in comparision with LHRS
// trying to see RHRRS RF-S2 time for the second time.
// Following exactely same steps from LHRS ie pad_7.c 
// Author Bishnu 
// Modify Itabashi 09/27/2018
// Pad# 7

double min(int i,int j);
double max(int i,int j);
void  hrs_tof_correct ()
{

  
 
  gStyle->SetOptStat(111111);
 TChain *T = new TChain("T");
 for(Int_t i = 100670; i<100675;i++) /// these run goes 100200 to 100685
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/bishnu/Rootfiles/tritium_%d*.root",i));

 }
 Int_t pmtl = 16; //PMT Left variable
 Int_t pmtr = 48;  // PMT Right variabl
 int ch; //channel
int tr; //triger conditon
 tr=1; //[0]: x correction [1]:y correction [2]:theata correction [3]: phi correction 
 ch=7; //select channel

 //----------- Comment out --------------------------------//
 cout<<"============= < Right arm S2 TOF analysis > ==================="<<endl;
 // cout<<"Channel number: "<<endl;
 // cin>>ch; 
 // cout<<"correction condition (x->0),(y->1),(theta->2),(phi->3):"<<endl;
 //cin>>tr;

 //------------------------------------------------------//

 int chmax=16;
 //====== Defined Function ================//
 double p0[chmax][4];
 double p1[chmax][4];
 TCanvas *c[chmax][4];
 TCut cut[chmax];
 TH1F *htof[chmax];
 TCanvas* cc[chmax];
 //=== X correction ====//
 TH2F *hx[chmax][4];
 TH2F *hx2[chmax][4];
 TH2F *hx3[chmax][4];
 TProfile *tpx[chmax][4];
 TProfile *tpx_it[chmax][4];
 TF1 *fx[chmax][4];
 TH2F *hx_c[chmax][4];
 TH1F *htof_xc[chmax][4];
 TH1F *htof_c[chmax];

 char *trx[4];
 trx[0]="R.tr.x[0]";
 trx[1]="R.tr.y[0]";
 trx[2]="R.tr.th[0]";
 trx[3]="R.tr.ph[0]";

 //===================== Cut condition===========================//

TCut cut_tof[chmax];
double tof_min[chmax],tof_max[chmax];

 tof_min[7]=0.224;
 tof_max[7]=0.228;


 cut_tof[7]=Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)>(%lf)*1.0e-6 && ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)<(%lf)*1.0e-6 ",pmtl+7,pmtr+7,tof_min[7],pmtl+7,pmtr+7,tof_max[7]);




 for(int i=0;i<chmax;i++){
 for(int j=0;j<4;j++){
   p0[i][j]=2.0e-9;
   p1[i][j]=2.0e-9;
}
}


 for(int i=0;i<chmax;i++){ 
if(i==ch ){
   for(int j=0;j<4;j++){
     //if(i==ch && j==tr ){     
     

       c[i][j] = new TCanvas(Form("c[%d][%d]",i,j),Form("c[%d][%d]",i,j));
   c[i][j]->Divide(2,2);

 cut[i]=(Form("DR.evtypebits>>7&1 &&  R.s2.t_pads[0]==%d    &&  R.s2.nthit==1",i));
 double xmin=min(i,j); 
 double xmax=max(i,j);
 //-------- TOF hist w/o calibration -----------//
 htof[i] = new TH1F(Form("htof[%d]",i),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6);
T->Project(Form("htof[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",pmtl+i,pmtr+i),cut[i]);

c[i][j]->cd(1);
 htof[i]->Draw();
  


//-------- TOF vs X hist --------------------------//
 hx[i][j] = new TH2F(Form("hx[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,xmin,xmax);


 T->Project(Form("hx[%d][%d]",i,j),Form("(%s):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0)",trx[j],pmtl+i,pmtr+i),cut[i] && cut_tof[i]);
c[i][j]->cd(2);
 hx[i][j]->Draw();


//---------TProfile  Correction for X ----------------------//



 hx2[i][j] = new TH2F(Form("hx2[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),200,xmin,xmax,300,0.2e-6,0.25e-6);
 T->Project(Form("hx2[%d][%d]",i,j),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0):(%s)",pmtl+i,pmtr+i,trx[j]),cut[i] && cut_tof[i]);

 c[i][j]->cd(3);
 hx2[i][j]->Draw();

 tpx[i][j] = new TProfile(Form("tpx[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 xposi ; ",i),200,xmin,xmax,0.2e-6,0.25e-6);
			  //,tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 
 T->Project(Form("tpx[%d][%d]",i,j),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0): %s",pmtl+i,pmtr+i,trx[j]),cut[i]);

  tpx[i][j]->Draw("same");


 
 //------ Fit TProfile Correction for X -------------------//

 fx[i][j]= new TF1(Form("fx[%d][%d]",i,j),"[0]+[1]*x",tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 fx[i][j]->SetParameter(0,p0[i][j]);
  fx[i][j]->SetParameter(1,p1[i][j]);
  tpx[i][j]->Fit(Form("fx[%d][%d]",i,j),"","");
 tpx[i][j]->GetYaxis()->SetRangeUser(tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 fx[i][j]->Draw("same"); 
 

 //----- Corrected for X -------------//  
 
 
 p0[i][j]=fx[i][j]->GetParameter(0);
 p1[i][j]=fx[i][j]->GetParameter(1);
 cout<<"Get parameter [0] is "<<p0[i][j]<<endl;
 cout<<"Get parameter [1] is "<<p1[i][j]<<endl;
 p1[i][j]=p1[i][j]*1.0e9; // [ns]



 /*

 //=================================//
//========= Iteration =============//
//===================================//

 tpx_it[i][j] = new TProfile(Form("tpx_it[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 xposi ; ",i),200,xmin,xmax,0.2e-6,0.25e-6);
 hx3[i][j] = new TH2F(Form("hx3[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),200,xmin,xmax,300,0.2e-6,0.25e-6);
   TF1* fx_it= new TF1("fx_it","[0]+[1]*x",0.2e-6,0.25e-6);
		       //tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
int s=0;
 while (s<3){

 T->Project(Form("hx3[%d][%d]",i,j),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - (RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0-(%s)*(1.0*(%lf*1.0e-9)): (%s)",pmtl+i,pmtr+i,trx[j],p1[i][j],trx[j]),cut[i]);

 c[i][j]->cd(4);
 hx3[i][j]->Draw();

 T->Project(Form("tpx_it[%d][%d]",i,j),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - (RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +(RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12)/2.0-(%s)*(1.0*(%lf*1.0e-9)): (%s)",pmtl+i,pmtr+i,trx[j],p1[i][j],trx[j]),cut[i]);


  
 tpx[i][j] = new TProfile(Form("tpx[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 xposi ; ",i),200,xmin,xmax,0.2e-6,0.25e-6);
			  //,tof_min[i]*1.0e-6,tof_max[i]*1.0e-6);
 
 T->Project(Form("tpx[%d][%d]",i,j),Form("((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0): %s",pmtl+i,pmtr+i,trx[j]),cut[i]);

 // -(%s)*(1.0*(%lf*1.0e-9)))",trx[j],pmtl+i,pmtr+i,trx[j],p1[i][j]),cut[i] &&cut_tof[i]);
  // tpx_it[i][j]->Fit(Form("fx_it",i,j),"","");

   tpx_it[i][j]->Fit("fx_it","","");
   //  c[i][j]->cd(4);
  tpx_it[i][j]->Draw("same");
 p0[i][j]=fx_it->GetParameter(0);
 p1[i][j]=fx_it->GetParameter(1);
 cout<<"Get parameter [0] is "<<p0[i][j]<<endl;
 cout<<"Get parameter [1] is "<<p1[i][j]<<endl;
  p1[i][j]=p1[i][j]*1.0e9; // [ns]
 cout<<"p1[i][j]: ("<<s<<")"<<p1[i][j]<<endl;
 s=s+1;
   } // end while 
*/


  c[i][j]->cd(4);
 hx_c[i][j] = new TH2F(Form("hx_c[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6,200,xmin,xmax);

 //p1[i][j]=8.0;

 T->Project(Form("hx_c[%d][%d]",i,j),Form("(%s):((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -(%s)*(1.0*(%lf*1.0e-9)))",trx[j],pmtl+i,pmtr+i,trx[j],p1[i][j]),cut[i] &&cut_tof[i]);
 hx_c[i][j]->Draw();


 
 //--------- TOF affter corrected X ----------//
 // p1[i][j]=20.23538;
  htof_xc[i][j] = new TH1F(Form("htof_xc[%d][%d]",i,j),Form("pad #%d RF -s2 time(s2 not corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6);

  T->Project(Form("htof_xc[%d][%d]",i,j),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -(%lf*1.0e-9)*%s)",pmtl+i,pmtr+i,p1[i][j],trx[j]),cut[i] && cut_tof[i]);


 c[i][j]->cd(1);
 htof_xc[i][j]->SetLineColor(2);
 htof_xc[i][j]->Draw("same");
 

 }

 cc[i] =new TCanvas(Form("cc[%d]",i),Form("cc[%d]",i));
 htof_c[i] = new TH1F(Form("htof_c[%d]",i),Form("pad #%d RF -s2 time(s2 all corrected) ;RF- s2 timein sec  ; ",i),300,0.2e-6,0.25e-6);
 T->Project(Form("htof_c[%d]",i),Form(" ((RTDC.F1FirstHit[15] -RTDC.F1FirstHit[14])*56.23e-12 - ((RTDC.F1FirstHit[15]-RTDC.F1FirstHit[%d])*56.23e-12 +((RTDC.F1FirstHit[46]-RTDC.F1FirstHit[%d])*56.23e-12))/2.0 -1.0/%lf*1.0e-9*%s-1.0*%lf*1.0e-9*%s-1.0/%lf*1.0e-9*%s-1.0/%lf*1.0e-9*%s)",pmtl+i,pmtr+i,p1[i][0],trx[0],p1[i][1],trx[1],p1[i][2],trx[2],p1[i][3],trx[3]),cut[i]&& cut_tof[i]);

 htof[i]->Draw();
 htof_c[i]->SetLineColor(kRed);
 htof_c[i]->Draw("same");


 }
 }

 


}

//================== X range function =================================//
double min(int i, int j){
  int chmax=30;
  double min[chmax][chmax];
  double mini;



  for(int k=0;k<chmax;k++){
    for(int l=0;l<chmax;l++){
      //defolt parameter
	min[k][l]=-1.0;
      

 }
  }

  min[7][0]=-0.25;
  min[7][1]=-0.05;
  min[7][2]=-0.05;
  min[7][3]=-0.05;
  for(int k=0;k<chmax;k++){
    for(int l=0;l<chmax;l++){
      if(i==k && j==l){
	mini=min[k][l];
      }   

 }
  }
  return mini;

}

double max(int i, int j){
  int chmax=30;
 double max[chmax][chmax];
  double maxi;
 
  

  for(int k=0;k<chmax;k++){
    for(int l=0;l<chmax;l++){
      //defolt parameter
	max[k][l]=1.0;
 }
  }



  max[7][0]=-0.02;
  max[7][1]=0.00;
  max[7][2]=0.0;
  max[7][3]=0.05;
  for(int k=0;k<chmax;k++){
    for(int l=0;l<chmax;l++){
      if(i==k && j==l){
	maxi=max[k][l];
      }

 }
  }

  return maxi;

}
