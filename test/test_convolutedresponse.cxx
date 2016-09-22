#include "WireCellSignal/FieldRsp.h"
#include "WireCellSignal/GenElecRsp.h"
#include "WireCellSignal/ConvolutedResponse.h"
#include "WireCellSignal/ElectronicsConfig.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include <iostream>

using namespace WireCellSignal;
using namespace std;

int main(int argc, char * argv[])
{
  
  ElectronicsConfig *EConfig = new ElectronicsConfig();
  ConvolutedResponse *cRsp = new ConvolutedResponse(EConfig);
  cRsp->OutputConvolutedResponse();     
  
  //int n = atoi(argv[1]);
  //ConvolutedResponse *cRsp = new ConvolutedResponse("/home/xiaoyueli/BNLIF/wire-cell/signal/convoluted_response.root");
  /*
  TCanvas c("c","",3200,1600);
  c.Divide(18,7);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	c.cd(18*k + i*6 + j +1);
	TH1F *htmp = cRsp->GetResponseFunction(i,j,k-3,n);
	htmp->Draw();
      }
    }
  }
  
  c.SaveAs(Form("%s_%d.pdf", "convoluted_response",n));  
  */
  return 0;
}
 
