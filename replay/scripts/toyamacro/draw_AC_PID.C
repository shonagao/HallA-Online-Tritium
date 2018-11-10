#include "Setting.cc"
void draw_AC_PID(){
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);


  //TFile *ifp = new TFile("pdf/test.root");
  //TFile *ifp = new TFile("pdf/Lambda_coin00.root");
  TFile *ifp = new TFile("pdf/Lambda_coin.root");
  Setting *set = new Setting();

  TLine *l_pi, *l_k, *l_p;
  l_pi = new TLine(0,0,0,0);
  l_k = new TLine(0,0,0,0);
  l_p = new TLine(0,0,0,0);
  set -> SetTLine(l_pi, 6, 1, 2);
  set -> SetTLine(l_k, 1, 1, 2);
  set -> SetTLine(l_p, 2, 1, 2);

  //time gate params
  double pi_min, pi_max, k_min, k_max, p_min, p_max;
  pi_min = -1.; pi_max = 1.;
  k_min  =  2.; k_max  = 4.4;
  p_min  =  9.; p_max  = 14.;

  TH2F *h2_s2ccoin_f1_a1_wa2cut, *h2_s2ccoin_f1_a2_wa1cut;
  TH1D *h_a1, *h_a2;
  TH1D *h_pi_a1, *h_pi_a2, *h_k_a1, *h_k_a2, *h_p_a1,*h_p_a2;

  h2_s2ccoin_f1_a1_wa2cut = (TH2F*)ifp->Get("h2_s2ccoin_f1_a1_wa2picut");
  h2_s2ccoin_f1_a2_wa1cut = (TH2F*)ifp->Get("h2_s2ccoin_f1_a2_wa1cut");


  double a1_min = -100.;
  double a1_max = 4000.;
  int pi_a1_min =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(pi_min);
  int pi_a1_max =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(pi_max);
  int k_a1_min  =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(k_min);
  int k_a1_max  =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(k_max);
  int p_a1_min  =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(p_min);
  int p_a1_max  =h2_s2ccoin_f1_a1_wa2cut ->GetXaxis()->FindBin(p_max);

  h_a1    = h2_s2ccoin_f1_a1_wa2cut -> ProjectionY("h_a1");
  h_pi_a1 = h2_s2ccoin_f1_a1_wa2cut -> ProjectionY("h_pi_a1",pi_a1_min, pi_a1_max);
  h_k_a1  = h2_s2ccoin_f1_a1_wa2cut -> ProjectionY("h_k_a1" ,k_a1_min , k_a1_max);
  h_p_a1  = h2_s2ccoin_f1_a1_wa2cut -> ProjectionY("h_p_a1" ,p_a1_min , p_a1_max);

  set -> SetTH1(h_a1     , "A1 ADC sum", "ADC [arb.]", "Counts", 1, 3000, 0);
  set -> SetTH1(h_pi_a1  , "#pi gated" , "", "", 1, 3001, 6);
  set -> SetTH1(h_k_a1   , "K gated"   , "", "", 1, 3001, 4);
  set -> SetTH1(h_p_a1   , "p gated"   , "", "", 1, 3001, 2);

  double a2_min = -100.;
  double a2_max = 20000.;
  int pi_a2_min =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(pi_min);
  int pi_a2_max =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(pi_max);
  int k_a2_min  =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(k_min);
  int k_a2_max  =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(k_max);
  int p_a2_min  =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(p_min);
  int p_a2_max  =h2_s2ccoin_f1_a2_wa1cut ->GetXaxis()->FindBin(p_max);

  h_a2    = h2_s2ccoin_f1_a2_wa1cut -> ProjectionY("h_a2");
  h_pi_a2 = h2_s2ccoin_f1_a2_wa1cut -> ProjectionY("h_pi_a2",pi_a2_min, pi_a2_max);
  h_k_a2  = h2_s2ccoin_f1_a2_wa1cut -> ProjectionY("h_k_a2" ,k_a2_min , k_a2_max);
  h_p_a2  = h2_s2ccoin_f1_a2_wa1cut -> ProjectionY("h_p_a2" ,p_a2_min , p_a2_max);

  set -> SetTH1(h_a2     , "A2 ADC sum", "ADC [arb.]", "Counts", 1, 3000, 0);
  set -> SetTH1(h_pi_a2  , "#pi gated" , "", "", 1, 3001, 6);
  set -> SetTH1(h_k_a2   , "K gated"   , "", "", 1, 3001, 4);
  set -> SetTH1(h_p_a2   , "p gated"   , "", "", 1, 3001, 2);

  TCanvas *c[4];
  for(int i=0;i<4;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,900);
  }



  c[0]->Clear();
  c[0]->cd(1);gPad->SetLogz(1); h2_s2ccoin_f1_a1_wa2cut->Draw("colz");
                                l_pi -> DrawLine(pi_min,a1_min ,pi_min,a1_max);
                                l_pi -> DrawLine(pi_max,a1_min ,pi_max,a1_max);
                                l_k  -> DrawLine(k_min ,a1_min ,k_min ,a1_max);
                                l_k  -> DrawLine(k_max ,a1_min ,k_max ,a1_max);
                                l_p  -> DrawLine(p_min ,a1_min ,p_min ,a1_max);
                                l_p  -> DrawLine(p_max ,a1_min ,p_max ,a1_max);

  c[1]->Clear();
  c[1]->Divide(1,4);
  c[1]->cd(1);  gPad->SetLogy(1);  h_a1    ->Draw("");
  c[1]->cd(2);  gPad->SetLogy(1);  h_pi_a1 ->Draw("");
  c[1]->cd(3);  gPad->SetLogy(1);  h_k_a1  ->Draw("");
  c[1]->cd(4);  gPad->SetLogy(1);  h_p_a1  ->Draw("");

  c[2]->Clear();
  c[2]->cd(1);gPad->SetLogz(1); h2_s2ccoin_f1_a2_wa1cut->Draw("colz");
                                l_pi -> DrawLine(pi_min,a2_min ,pi_min,a2_max);
                                l_pi -> DrawLine(pi_max,a2_min ,pi_max,a2_max);
                                l_k  -> DrawLine(k_min ,a2_min ,k_min ,a2_max);
                                l_k  -> DrawLine(k_max ,a2_min ,k_max ,a2_max);
                                l_p  -> DrawLine(p_min ,a2_min ,p_min ,a2_max);
                                l_p  -> DrawLine(p_max ,a2_min ,p_max ,a2_max);

  c[3]->Clear();
  c[3]->Divide(1,4);
  c[3]->cd(1);  gPad->SetLogy(1);  h_a2    ->Draw("");
  c[3]->cd(2);  gPad->SetLogy(1);  h_pi_a2 ->Draw("");
  c[3]->cd(3);  gPad->SetLogy(1);  h_k_a2  ->Draw("");
  c[3]->cd(4);  gPad->SetLogy(1);  h_p_a2  ->Draw("");


  c[0]->Print("pdf/AC_PIC_Lambda.pdf[");
  c[0]->Print("pdf/AC_PIC_Lambda.pdf");
  c[1]->Print("pdf/AC_PIC_Lambda.pdf");
  c[2]->Print("pdf/AC_PIC_Lambda.pdf");
  c[3]->Print("pdf/AC_PIC_Lambda.pdf");
  c[3]->Print("pdf/AC_PIC_Lambda.pdf]");

}
