// 10/4/2018 Author Itabshi
// time walk correction R-HRS S2

void hrs_tof_wt_correction(){


 
  gStyle->SetOptStat(000000);

  //-------- TTree data input ---------------//
 TChain *T = new TChain("T");
 for(Int_t i = 93500; i<93501;i++) /// these run goes 93495 to 93504
 {
 T->Add(Form("/adaqfs/home/a-onl/tritium_work/itabashi/HallA-Online-Tritium/replay/Rootfiles/pass1/tritium_%d*.root",i));
 }

 //========== Defined Parameters =======================//
 int chmax=16; // channel of S2 PMT
 // chmax=1; //test 
 int ch=8;

 TH2F* htof_adc[chmax][2];
 TCanvas* c[chmax];
 TCanvas* c0[chmax];
 TProfile* tp[chmax][2];
 TF1* ffit[chmax][2];
 TH2F* htdc_adc[chmax][2];
 double tdc=0.5;//TDC converse ch->sec [ns/ch]
 double min_tof,max_tof,min_adc,max_adc;
 int bin_tof,bin_adc;
 //min_tof=-2.0e-8;
 //max_tof=2.0e-8;
 min_tof=2.0e-7;
 max_tof=2.2e-7;
 bin_tof=8000;
 min_adc=500.;
 max_adc=5000.;
 bin_adc=10000.;
 double p0[chmax][2],p1[chmax][2];


 for(int i=0;i<chmax;i++){
   

     if(i==ch){

   for(int j=0;j<2;j++){
     htof_adc[i][j]=new TH2F(Form("htof_adc[%d][%d]",i,j),Form("TOF vs S2 CH %d,[%d] PMT adc correlation ",i,j),bin_adc,min_adc,max_adc,bin_tof,min_tof,max_tof);

     //--------- TDC vs ADC ---------------//
  htdc_adc[i][j]=new TH2F(Form("htof_adc[%d][%d]",i,j),Form("TOF vs S2 CH %d,[%d] PMT adc correlation ",i,j),bin_adc,min_adc,max_adc,bin_tof,min_tof,max_tof);
      }


   //     c[i]=new TCanvas(Form("c[%d]",i),Form("c[%d]",i));
     c0[i]=new TCanvas(Form("c0[%d]",i),Form("c0[%d]",i));

     //------ Project hist -----------------------//
     //     T->Project(Form("htof_adc[%d][0]",i),Form("((R.s2.lt_c[%d]+R.s2.rt_c[%d])/2.0-(R.s0.lt_c+R.s0.rt_c)/2.0):R.s2.ra_c[%d]",i,i,i));//Right PMT
     //   T->Project(Form("htof_adc[%d][1]",i),Form("((R.s2.lt_c[%d]+R.s2.rt_c[%d])/2.0-(R.s0.lt_c+R.s0.rt_c)/2.0):R.s2.la_c[%d]",i,i,i));//Left PMT

   T->Project(Form("htdc_adc[%d][0]",i),Form("R.s2.rt_c[%d]:R.s2.ra_c[%d]",i,i));//Right PMT
     T->Project(Form("htdc_adc[%d][1]",i),Form("R.s2.lt_c[%d]:R.s2.la_c[%d]",i,i));//Left PMT



     //-------------- TProfile Function-----------------//
     tp[i][0]=new TProfile(Form("tp[%d][0]",i),Form("TProfile S2 ch %d R-PMT",i),bin_adc/20,min_adc,max_adc,min_tof,max_tof);
   tp[i][1]=new TProfile(Form("tp[%d][1]",i),Form("TProfile S2 ch %d L-PMT",i),bin_adc/10,min_adc,max_adc,min_tof,max_tof);



  //------- TDC Prjection ----------//
    T->Project(Form("tp[%d][0]",i),Form("R.s2.rt_c[%d]:R.s2.ra_c[%d]",i,i));//Right PMT
    T->Project(Form("tp[%d][1]",i),Form("R.s2.lt_c[%d]:R.s2.la_c[%d]",i,i));//Left PMT


   //------- TOF Prjection ----------//
   //  T->Project(Form("tp[%d][0]",i),Form("((R.s2.lt_c[%d]+R.s2.rt_c[%d])/2.0-(R.s0.lt_c+R.s0.rt_c)/2.0):R.s2.ra_c[%d]",i,i,i));//Right PMT
   //  T->Project(Form("tp[%d][1]",i),Form("((R.s2.lt_c[%d]+R.s2.rt_c[%d])/2.0-(R.s0.lt_c+R.s0.rt_c)/2.0):R.s2.la_c[%d]",i,i,i));//Left PMT

     //------------ Fit Function -------------------//
     ffit[i][0]=new TF1(Form("ffit[%d][0]",i),"[0]*(1./sqrt(x)-1./sqrt([1]))",min_adc,max_adc);
     //  ffit[i][0]=new TF1(Form("ffit[%d][0]",i),"[0]*(1./sqrt(x)-[1])",100.,1000.);
     ffit[i][1]=new TF1(Form("ffit[%d][1]",i),"[0]*(pow(x,-1/2)-pow([1],-1/2))",min_adc,max_adc);
     // ffit[i][0]->SetParameter(1,);
     tp[i][0]->Fit(Form("ffit[%d][0]",i),"","");
     tp[i][1]->Fit(Form("ffit[%d][1]",i),"","");
     p0[i][0]=ffit[i][0]->GetParameter(0);
     p0[i][1]=ffit[i][1]->GetParameter(0);
     p1[i][0]=ffit[i][0]->GetParameter(1);
     p1[i][1]=ffit[i][1]->GetParameter(1);

    //----------------- Comment out -----------------//
		      cout<<"======== S2 R-PMT ch"<<i<<"parameters ========"<<endl;
		      cout<<"p0 :"<<p0[i][0]<<"  p1 :"<<p1[i][0]<<endl;    
                      cout<<"======== S2 L-PMT ch"<<i<<"parameters ========"<<endl;
		      cout<<"p0 :"<<p0[i][1]<<"  p1 :"<<p1[i][1]<<endl;    


     //-----------TCanvas -----------------------//
		    

      //----- TOF Canvas ---------//
		      /*
			 c[i]->Divide(1,2);
			 c[i]->cd(1);
			 htof_adc[i][0]->Draw("colz");
			 tp[i][0]->Draw("same");
                         c[i]->cd(2);
			 htof_adc[i][1]->Draw("colz");
			  tp[i][1]->Draw("same");
		      */
		      //----- ADC vs TDC ----------//
			 c0[i]->Divide(1,2);
			 c0[i]->cd(1);
			 htdc_adc[i][0]->Draw("colz");
			 tp[i][0]->Draw("same");
                         c0[i]->cd(2);
			 htdc_adc[i][1]->Draw("colz");
			 tp[i][1]->Draw("same");

			 }
 }





}
