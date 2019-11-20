#include "WCPSignal/GenNoise.h"
#include "WCPSignal/ElectronicsConfig.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace WCPSignal;
using namespace std;

int main(int argc, char * argv[])
{
  ElectronicsConfig *EConfig = new ElectronicsConfig();
  //EConfig->SetShapingTime(1);
  //EConfig->SetGain(7.8);
  GenElecRsp *ns = new GenElecRsp(EConfig);
  TCanvas *c = new TCanvas();
  TF1 *shp = (TF1*)ns->GetShapingFunction();
  shp->SetLineColor(kRed);
  c->cd();
  shp->Draw();
  c->SaveAs("default_eRsp.pdf");
}
