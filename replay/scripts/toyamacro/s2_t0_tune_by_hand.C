#include "ParamMan.cc"
#include "define.h"

//This macro is useful when you change peak offset of cointime.
//If you want to move peak position +3ns, you can set offsetL=-3.
//Then peak pos will move to positive side
void s2_t0_tune_by_hand(){

  ParamMan *param = new ParamMan("param/f1_tuned_Lambda_twc.param");
  //ParamMan *param = new ParamMan("param/f1_tuned_Lambda.param");
  param->SetVal();

  double offsetL =-3.;
  int LR;

  for(int i= 0;i<16;i++){
    LR = 0; 
    param->SetTimeTune(CID_F1S2, i, LR, 0, offsetL);
    param->SetTimeTune(CID_F1S2, i, LR, 1, offsetL);
  }


  param->WriteToFile("param/tmp.param");
}
