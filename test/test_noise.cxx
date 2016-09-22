#include "WireCellSignal/GenNoise.h"
#include "WireCellSignal/ElectronicsConfig.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace WireCellSignal;
using namespace std;

int main(int argc, char * argv[])
{
  ElectronicsConfig *EConfig = new ElectronicsConfig();
  //EConfig->SetNTDC(1000);
  EConfig->SetShapingTime(2);
  EConfig->SetGain(7.8);
  GenNoise *ns = new GenNoise(*EConfig);
  std::cout<<"ENC: "<<ns->NoiseRMS()<<std::endl;
  std::cout<<"baseline is "<<ns->GetBaseline()<<std::endl;
  TH1 *f = (TH1*)ns->NoiseInTime();
  TCanvas *c = new TCanvas();
  ns->PrintNoiseInTime(c);
  //ns->PrintNoiseInFrequency(c);
  return 0;
}
