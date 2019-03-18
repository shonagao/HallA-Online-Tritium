/*
  "correction.hh"
  
  Toshi Gogami , 15Oct2014
*/


#include <TMath.h>

int filled_event = 0;
//const int max_event = 10000;
const int max_event = 30000;
//double a1[max_event];
//double a2[max_event];
//double a3[max_event];
//double a4[max_event];
int a1[max_event];
int a2[max_event];
int a3[max_event];
int a4[max_event];
double t1[max_event];
double t2[max_event];
double t3[max_event];
double t4[max_event];
double par[100]; // parameters
Int_t allparam = 0;
double w1, w2, w3, w4, wtof; // weights of chi-square

extern double calcTime(double,double,double,int,double);
extern void fcn(Int_t &nPar, Double_t* /*grad*/, Double_t &fval, Double_t* param, Int_t /*iflag*/);
extern double tune(double* pa, int j);
