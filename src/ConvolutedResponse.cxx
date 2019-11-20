#include "WCPSignal/ConvolutedResponse.h"

WCPSignal::ConvolutedResponse::ConvolutedResponse(ElectronicsConfig *conf, const std::string garfield_file, const float overall_time_offset, const float uv_time_offset, const float uw_time_offset, const std::string outfilename)
  : outfile(outfilename)
{
  fRsp = new FieldRsp(conf, garfield_file.c_str(), 0, overall_time_offset, uv_time_offset, uw_time_offset);
  eRsp = new GenElecRsp(conf);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	hRsp[i][j][k] = (TH1*)(fRsp->GetFieldResponseHist(i,j,k)->Clone(Form("hRsp_%d_%d_%d",i,j,k)));
	hRsp[i][j][k]->Reset();
	hRspMag[i][j][k] = (TH1*)(fRsp->GetFieldResponseInFreq(i,j,k)->Clone(Form("hRsp_mag_%d_%d_%d",i,j,k)));
	hRspPh[i][j][k] = (TH1*)(fRsp->GetFieldResponseInPhase(i,j,k)->Clone(Form("hRsp_ph_%d_%d_%d",i,j,k)));
      }
    }
  }
  BuildConvolutedResponse();
  h = (TH1F*)hRsp[0][0][0]->Clone();
  h->Reset();
}

WCPSignal::ConvolutedResponse::ConvolutedResponse(const std::string infilename)
{
  infile = new TFile(infilename.c_str(), "read");
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	hRsp[i][j][k] = (TH1*)(infile->Get(Form("hRsp_%d_%d_%d",i,j,k))->Clone());
	if (!hRsp[i][j][k]) {
	  std::cerr << "cannot find convoluted response function for plane "
		    << i << ", point " << j << ", wire " << k << std::endl;
	  exit(-1);
	}
	hRspMag[i][j][k] = hRsp[i][j][k]->FFT(0,"MAG");
	hRspPh[i][j][k] = hRsp[i][j][k]->FFT(0,"PH");
      }
    }
  }
  h = (TH1F*)hRsp[0][0][0]->Clone();
  h->Reset();
}

WCPSignal::ConvolutedResponse::ConvolutedResponse(FieldRsp *fR, GenElecRsp *eR, const std::string outfilename)
  : fRsp(fR)
  , eRsp(eR)
  , outfile(outfilename)
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	hRsp[i][j][k] = (TH1*)(fRsp->GetFieldResponseHist(i,j,k)->Clone(Form("hRsp_%d_%d_%d",i,j,k)));
	hRsp[i][j][k]->Reset();
	hRspMag[i][j][k] = (TH1*)(fRsp->GetFieldResponseInFreq(i,j,k)->Clone(Form("hRsp_mag_%d_%d_%d",i,j,k)));
	hRspPh[i][j][k] = (TH1*)(fRsp->GetFieldResponseInPhase(i,j,k)->Clone(Form("hRsp_ph_%d_%d_%d",i,j,k)));
      }
    }
  }
  BuildConvolutedResponse();
  h = (TH1F*)hRsp[0][0][0]->Clone();
  h->Reset();
  //h->Rebin(5);
}

WCPSignal::ConvolutedResponse::~ConvolutedResponse()
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	if (hRsp[i][j][k]) hRsp[i][j][k]->Delete();
	if (hRspMag[i][j][k]) hRspMag[i][j][k]->Delete();
	if (hRspPh[i][j][k]) hRspPh[i][j][k]->Delete();
      }
    }
  }
  if (infile) infile->Close();
  if (eRsp) delete eRsp;
  if (fRsp) delete fRsp;
}

TH1F* WCPSignal::ConvolutedResponse::GetResponseFunction(int plane, int index1, int index2, int Tshift, double charge)
{
  //*h = (TH1F*)hRsp[plane][index1][index2+3];//->Clone();
  std::cout<<"plane = "<<plane<<", index1 = "<<index1<<", index2 = "<<index2<<std::endl;
  h->SetName(Form("h_%d_%d_%d",plane,index1,index2+3));
  h->Reset();
  if (plane < 0 || plane > 2 ||
      index1 < 0 || index1 > 5 ||
      index2 < -3 || index2 > 3) {
    std::cerr << "plane = "<<plane
	      <<" index1 = "<<index1
	      <<" index2 = "<<index2
      << "... wrong input to function" << std::endl;
    //exit(-1);
    return h;
  }  
  for (int i = 0; i < hRsp[plane][index1][index2+3]->GetNbinsX(); ++i) {
    if (i<Tshift) {
      h->SetBinContent(i+1, 0);
    } else {
      h->SetBinContent(i+1, hRsp[plane][index1][index2+3]->GetBinContent(i+1-Tshift)); 
    }
  }
  h->Scale(charge);
  return h;  
}

void WCPSignal::ConvolutedResponse::BuildConvolutedResponse()
{
  int NTDC = hRsp[0][0][0]->GetNbinsX();
  double rho, phi;
  TVirtualFFT *ifft;
  double *re, *im;
  re = new double[NTDC];
  im = new double[NTDC];
  
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	for (int b = 1; b < NTDC+1; ++b) {
	  rho = hRspMag[i][j][k]->GetBinContent(b) * eRsp->GetShapingFInFreq()->GetBinContent(b);
	  phi = hRspPh[i][j][k]->GetBinContent(b) + eRsp->GetShapingFInPhase()->GetBinContent(b);
	  hRspMag[i][j][k]->SetBinContent(b, rho);
	  hRspPh[i][j][k]->SetBinContent(b, phi);
	  re[b-1] = rho * TMath::Cos(phi) / NTDC;
	  im[b-1] = rho * TMath::Sin(phi) / NTDC;
	}
	ifft = TVirtualFFT::FFT(1, &NTDC, "C2R M K");
	ifft->SetPointsComplex(re,im);
	ifft->Transform();
	hRsp[i][j][k] = TH1::TransformHisto(ifft, 0, "Re");
	hRsp[i][j][k]->SetName(Form("hRsp_%d_%d_%d",i,j,k));
	hRsp[i][j][k]->Rebin(5);
	hRsp[i][j][k]->Scale(0.2);
      }
    }
  }
  delete re;
  delete im;
}
/*
TH1* WCPSignal::ConvolutedResponse::GetResponseMag(int plane, int point, int wire)
{  
  return hRspMag[plane][abs(point-4)][wire];
}

TH1* WCPSignal::ConvolutedResponse::GetResponsePh(int plane, int point, int wire)
{
  return hRspPh[plane][abs(point-4)][wire];
}
*/
void WCPSignal::ConvolutedResponse::OutputConvolutedResponse()
{
  TFile f(outfile.c_str(), "recreate");
  f.cd();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	hRsp[i][j][k]->Write();
      }
    }
  }
  f.Close();
}
