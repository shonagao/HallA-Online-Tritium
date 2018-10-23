/*
  "ParamMan.hh"
  
  Toshi Gogami , 15Oct2014
*/

#ifndef ParamMan_h
#define ParamMan_h 1

#include <iostream>
using namespace std;

class ParamMan{
  
public:
  ParamMan();
  ~ParamMan();
  
private:
  double par1[3];
  double par2[3];
  double par3[3];
  double par4[3];
  
public:
  void SetParam1(double p1,double p2, double p3){
    par1[0] = p1;
    par1[1] = p2;
    par1[2] = p3;
  };
  void SetParam2(double p1,double p2, double p3){
    par2[0] = p1;
    par2[1] = p2;
    par2[2] = p3;
  };
  void SetParam3(double p1,double p2, double p3){
    par3[0] = p1;
    par3[1] = p2;
    par3[2] = p3;
  };
  void SetParam4(double p1,double p2, double p3){
    par4[0] = p1;
    par4[1] = p2;
    par4[2] = p3;
  };
  
  double* GetParam1(){
    return par1;
  };
  double* GetParam2(){
    return par2;
  };
  double* GetParam3(){
    return par3;
  };
  double* GetParam4(){
    return par4;
  };
  
  void CheckParam(){
    cout << "CH1: " << par1[0] << " " << par1[1] << " " << par1[2] << endl;
    cout << "CH2: " << par2[0] << " " << par2[1] << " " << par2[2] << endl;
    cout << "CH3: " << par3[0] << " " << par3[1] << " " << par2[2] << endl;
    cout << "CH4: " << par4[0] << " " << par4[1] << " " << par3[2] << endl;
  }
  
  
};

#endif
