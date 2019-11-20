#include "WCPSignal/FieldRsp.h"
#include "WCPSignal/ElectronicsConfig.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace WCPSignal;
using namespace std;

int main(int argc, char * argv[])
{
  ElectronicsConfig *EConfig = new ElectronicsConfig();
  TFile *outfile = new TFile("field_response.root","recreate");
  FieldRsp *fieldR = new FieldRsp(EConfig, "/home/xiaoyueli/BNLIF/wire-cell/signal/dune.root", outfile);
  //fieldR->SetOutfile(outfile);
  fieldR->DrawFieldResponse();
  //fieldR->DrawFieldResponseFFT();
}
