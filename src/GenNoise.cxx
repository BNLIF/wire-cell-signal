#include "WireCellSignal/GenNoise.h"

double PoissonReal(const Double_t *k,  const Double_t *lambda) {
  return TMath::Exp(k[0]*TMath::Log(lambda[0])-lambda[0]) / TMath::Gamma(k[0]+1.);
}

WireCellSignal::GenNoise::GenNoise(ElectronicsConfig& con, double rand)
  : config(&con)
{
  Init(config->ShapingTime());
  eRspD = new GenElecRsp(config, 14, 2);
  hERspD = eRspD->GetShapingFunctionHist();
  hERspD->Rebin(5);
  hERspD->Scale(0.2);
  TVirtualFFT::SetTransform(0);
  eRspMagD = hERspD->FFT(0, "MAG");
  eRspPhD = hERspD->FFT(0, "PH");
  eRspMagD->SetName("defaultMag");
  eRspPhD->SetName("defaultPh");
  //eRspMagD = (TH1*)eRspD->GetShapingFInFreq()->Clone("defaultMag");
  //eRspPhD = (TH1*)eRspD->GetShapingFInPhase()->Clone("defaultPh");

  eRsp = new GenElecRsp(config, config->Gain(), config->ShapingTime());
  hERsp= eRsp->GetShapingFunctionHist();
  hERsp->Rebin(5);
  hERsp->Scale(0.2);
  eRspMag = hERsp->FFT(0, "MAG");
  eRspPh = hERsp->FFT(0, "PH");
  //eRspMag = eRsp->GetShapingFInFreq();
  //eRspPh = eRsp->GetShapingFInPhase();

  hFreq = new TH1F("hFreq", "", config->NTDC(), 0, config->NTDC());
  
  re = new double[config->NTDC()];
  im = new double[config->NTDC()];
  fb = 0;
  PinkNoise(rand);
  baseline = -1;
}

void WireCellSignal::GenNoise::Init(double shapingtime)
{
  double MaxPoissArg = 100.;
  MyPoisson = new TF1("MyPoisson", PoissonReal, 0., MaxPoissArg, 1);
  MyPoisson->SetParameter(0, 3.3708);
  digitNoise = 4.27132e+01;
  f1 = new TF1("f","([0]+[1]*x*9592./2.)*exp(-[2]*pow(x*9592./2.,[3]))"); // x in MHz
  double fitpar[4]={6.22750e+02,-2.53535e-01,8.07578e-05,1.35510e+00};// how to get these numbers?
  f1->SetParameters(fitpar);  
}

WireCellSignal::GenNoise::~GenNoise()
{
  delete f1;
  delete MyPoisson;
  delete hFreq;
  delete re;
  delete im;
}

void WireCellSignal::GenNoise::PinkNoise(double rand)
{
  int NTDC = config->NTDC();
  hFreq->Reset();
  for (int i = 0; i < NTDC; ++i) {
    double freq;
    if ( i < NTDC/2. ) freq = config->DigitFreq()*i/NTDC;
    else freq = config->DigitFreq() * (1.0 - (double)i/NTDC);
    double rho = f1->Eval(freq) * eRspMag->GetBinContent(i+1) / eRspMagD->GetBinContent(i+1) + digitNoise;
    if (!i) rho = 0;
    rho *= std::sqrt(NTDC/9592.) * MyPoisson->GetRandom()/3.3708;
    hFreq->SetBinContent(i+1, rho);
    double phi = gRandom->Uniform(0, 1) * 2. * TMath::Pi();
    phi += eRspPh->GetBinContent(i+1) - eRspPhD->GetBinContent(i+1);
    re[i] = rho * std::cos(phi) / NTDC;
    im[i] = rho * std::sin(phi) / NTDC;
  }
  TVirtualFFT *ifft = TVirtualFFT::FFT(1, &NTDC, "C2R M K");
  ifft->SetPointsComplex(re, im);
  ifft->Transform();
  fb = TH1::TransformHisto(ifft, fb, "Re");
  fb->SetName("noiseInTime");
}

double WireCellSignal::GenNoise::NoiseRMS()
{
  double rms = 0;
  for (int i = 0; i < config->NTDC(); ++i) {
    rms += std::pow(fb->GetBinContent(i+1), 2);
  }
  rms = std::sqrt(rms/config->NTDC());
  return rms;
}

void WireCellSignal::GenNoise::PrintNoiseInTime(TCanvas *c)
{
  if (c) c->cd();
  fb->Draw();
  c->SaveAs("noiseT.pdf");
}

void WireCellSignal::GenNoise::PrintNoiseInFrequency(TCanvas *c)
{
  if (c) c->cd();
  hFreq->Draw();
  c->SaveAs("noiseF.pdf");
}

TH1* WireCellSignal::GenNoise::NoiseInTime()
{
  fb->Reset();
  PinkNoise();
  return fb;
}

double WireCellSignal::GenNoise::GetBaseline()
{
  if (baseline>0) return baseline;
  else {
    baseline = 0;
    for(int i = 0; i < 1000; ++i) {
      fb->Reset();
      PinkNoise(gRandom->Uniform(0, 1));
      double n = NoiseRMS();
      baseline += n;      
    }
    baseline /= 1000;
    //    std::cout << "Noise RMS baseline is " << baseline << std::endl;
  }
  return baseline;
}
