#include "WCPSignal/FieldRsp.h"
#include "TString.h"

WCPSignal::FieldRsp::FieldRsp(ElectronicsConfig *con, const float overall_time_offset, const float uv_time_offset, const float uw_time_offset)
  : config(con)
{
  simu2D = false;
  outfile = 0;
  InitFieldResponse(overall_time_offset, uv_time_offset, uw_time_offset);
  FFTFieldResponse();
}

WCPSignal::FieldRsp::FieldRsp(ElectronicsConfig *con, const std::string infilename, TFile *f, const float overall_time_offset, const float uv_time_offset, const float uw_time_offset)
  : config(con)
  , outfile(f)
{
  simu2D = true;
  InitFieldResponse(infilename.c_str(), overall_time_offset, uv_time_offset, uw_time_offset);
  FFTFieldResponse();
}

WCPSignal::FieldRsp::~FieldRsp()
{
  if(!simu2D) {
    for (int i = 0; i < 3; ++i) {
      if (fieldResp[i]) delete fieldResp[i];
      if (gFieldResp[i]) delete gFieldResp[i];
      if (hFieldRespInFreq[i]) delete hFieldRespInFreq[i];
      if (hFieldRespInPhase[i]) delete hFieldRespInPhase[i];
    }
  } else {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  //if (i==2 && k!=3) continue;
	  if (hFieldResp2D[i][j][k]) delete hFieldResp2D[i][j][k];
	  if (hFieldRespInFreq2D[i][j][k]) delete hFieldRespInFreq2D[i][j][k];
	  if (hFieldRespInPhase2D[i][j][k]) delete hFieldRespInPhase2D[i][j][k];
	}
      }
    }
  }
  if(outfile) outfile->Close();
}

void WCPSignal::FieldRsp::InitFieldResponse(const float overall_time_offset, const float uv_time_offset, const float uw_time_offset)
{
#include "../../2dtoy/src/data.txt" // not ideal. shall be changed later
  gFieldResp[0] = new TGraph(5000, xu, yu);
  gFieldResp[1] = new TGraph(5000, xv, yv);
  gFieldResp[2] = new TGraph(5000, xw, yw);
  int nbins = config->NTDC() * 5;
  fieldResp[0] = new TH1F("fieldR_U", "fieldR_U", nbins, 0, nbins);
  fieldResp[1] = new TH1F("fieldR_V", "fieldR_V", nbins, 0, nbins);
  fieldResp[2] = new TH1F("fieldR_W", "fieldR_W", nbins, 0, nbins);
  for (int i = 0; i < nbins; ++i) {
    double time = fieldResp[0]->GetBinCenter(i+1) - 50. - overall_time_offset; // time in us
    fieldResp[0]->SetBinContent(i+1, gFieldResp[0]->Eval(time));
    fieldResp[1]->SetBinContent(i+1, gFieldResp[1]->Eval(time - uv_time_offset));
    fieldResp[2]->SetBinContent(i+1, gFieldResp[2]->Eval(time - uw_time_offset));
  }
  double scale = fieldResp[2]->Integral();
  fieldResp[0]->Scale(1./scale);
  fieldResp[1]->Scale(1./scale);
  fieldResp[2]->Scale(1./scale);
  if(outfile) {
    outfile->cd();
    gFieldResp[0]->Write();
    gFieldResp[1]->Write();
    gFieldResp[2]->Write();
    fieldResp[0]->Write();
    fieldResp[1]->Write();
    fieldResp[2]->Write();
  }
}

void WCPSignal::FieldRsp::InitFieldResponse(const std::string infilename, const float overall_time_offset, const float uv_time_offset, const float uw_time_offset)
{
  TFile infile(infilename.c_str(), "read");
  for (int j = 0; j < 6; ++j) {
    for (int k = 0; k < 7; ++k) {
      gFieldResp2D[0][j][k] = (TGraph*)infile.Get(Form("guc_%d",j*7+k+1));
      gFieldResp2D[1][j][k] = (TGraph*)infile.Get(Form("gvc_%d",j*7+k+1));
      if (k!=3) gFieldResp2D[2][j][k] = (TGraph*)infile.Get(Form("gwc_%d",j*7+k+1));
      else gFieldResp2D[2][j][k] = (TGraph*)infile.Get(Form("gwi_%d",j*7+k+1));
    }
    //gFieldResp2D[2][j][3] = (TGraph*)infile.Get(Form("gwi_%d",j*7+4));
  }
  infile.Close();
  
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 6; ++j) {
      for (int k = 0; k < 7; ++k) {
	hFieldResp2D[i][j][k] = new TH1F(Form("fR_%d_%d_%d",i,j,k), Form("fR_%d_%d_%d",i,j,k), config->NTDC()*5, overall_time_offset, config->NTDC()/config->DigitFreq()+overall_time_offset);
      }
    }
  }
  for(int k = 0; k < config->NTDC()*5; ++k) {
    double time = hFieldResp2D[0][0][0]->GetBinCenter(k+1) - overall_time_offset;
    for(int i = 0; i < 6; ++i){
      for (int j = 0; j < 7; ++j) {
	hFieldResp2D[0][i][j]->SetBinContent(k+1, gFieldResp2D[0][i][j]->Eval(time)); // U 
	hFieldResp2D[1][i][j]->SetBinContent(k+1, gFieldResp2D[1][i][j]->Eval(time-uv_time_offset)); // V
	/*if (j==3)*/
	hFieldResp2D[2][i][j]->SetBinContent(k+1, gFieldResp2D[2][0][j]->Eval(time-uw_time_offset)); // Y
	//else hFieldResp2D[2][i][j]->SetBinContent(k+1, 0); 
      }
    }
  }

  for(int i = 0; i < 6; ++i) {
    double scale = TMath::Abs(hFieldResp2D[2][i][3]->Integral());
    //hFieldResp2D[2][i][3]->Scale(1./scale);
    for (int k = 0; k < 7; ++k) {
      hFieldResp2D[0][i][k]->Scale(1./scale);
      hFieldResp2D[1][i][k]->Scale(1./scale);
      hFieldResp2D[2][i][k]->Scale(1./scale);
    }
  }

  if(outfile) {
    outfile->cd();
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  hFieldResp2D[i][j][k]->Write();
	}
      }
    }
  }
}

void WCPSignal::FieldRsp::DrawFieldResponse(const std::string canvas_name)
{
  TCanvas c(canvas_name.c_str(),"",32000,16000);
  if(!simu2D) {
    c.Divide(1,3);
    c.cd(1);
    fieldResp[0]->Draw();
    c.cd(2);
    fieldResp[1]->Draw();
    c.cd(3);
    fieldResp[2]->Draw();
  } else {
    c.Divide(18,7);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  //if (i==2 && k!=3) continue;
	  c.cd(18*k + i*6 + j +1); gFieldResp2D[i][j][k]->Draw("AC");
	}
      }
    }
    c.SaveAs(Form("%s_graph.pdf", canvas_name.c_str()));
    c.Clear(); c.Divide(18,7);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  c.cd(18*k + i*6 + j +1); hFieldResp2D[i][j][k]->Draw();
	}
      }
    }
  }
  c.SaveAs(Form("%s.pdf", canvas_name.c_str()));
}

void WCPSignal::FieldRsp::FFTFieldResponse()
{
  TVirtualFFT::SetTransform(0);
  TVirtualFFT *ifft;
  if (!simu2D) {
    for(int i = 0; i < 3; ++i) {
      hFieldRespInFreq[i] = fieldResp[i]->FFT(0, "MAG");
      //hFieldRespInFreq[i]->Scale(1./TMath::Sqrt((double)hFieldRespInFreq[i]->GetNbinsX()));
      hFieldRespInFreq[i]->SetName(Form("field_response_fft_mag_%d",i));
      hFieldRespInFreq[i]->SetTitle(Form("field_response_fft_mag_%d",i)); 
      hFieldRespInPhase[i] = fieldResp[i]->FFT(0, "PH");
      hFieldRespInPhase[i]->SetName(Form("field_response_fft_ph_%d",i));
      hFieldRespInPhase[i]->SetTitle(Form("field_response_fft_ph_%d",i));
      if(outfile) {
	outfile->cd();
	hFieldRespInFreq[i]->Write();
	hFieldRespInPhase[i]->Write();
      }
    }
  }
  else {
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  //if(i==2 && k!=3) continue;
	  hFieldRespInFreq2D[i][j][k] = hFieldResp2D[i][j][k]->FFT(0, "MAG");
	//hFieldRespInFreq2D[i][j]->Scale(1./TMath::Sqrt((double)hFieldRespInFreq2D[i][j]->GetNbinsX()));
	  hFieldRespInFreq2D[i][j][k]->SetName(Form("field_response_fft_mag_%d_%d_%d",i,j,k));
	  hFieldRespInFreq2D[i][j][k]->SetTitle(Form("field_response_fft_mag_%d_%d_%d",i,j,k)); 
	  hFieldRespInPhase2D[i][j][k] = hFieldResp2D[i][j][k]->FFT(0, "PH");
	  hFieldRespInPhase2D[i][j][k]->SetName(Form("field_response_fft_ph_%d_%d_%d",i,j,k));
	  hFieldRespInPhase2D[i][j][k]->SetTitle(Form("field_response_fft_ph_%d_%d_%d",i,j,k));
	  if(outfile) {
	    outfile->cd();
	    hFieldRespInFreq2D[i][j][k]->Write();
	    hFieldRespInPhase2D[i][j][k]->Write();
	  }
	}
      }
    }
  }
}

void WCPSignal::FieldRsp::DrawFieldResponseFFT(const std::string canvas_name)
{
  TCanvas c(canvas_name.c_str(),"",32000,32000);
  if(!simu2D) {
    c.Divide(2,3);
    for(int i = 0; i < 3; ++i) {
      c.cd(i*2+1);
      hFieldRespInFreq[i]->Draw();
      c.cd(i*2+2);
      hFieldRespInPhase[i]->Draw();
    }
  } else {
    c.Divide(13,14);
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 6; ++j) {
	for (int k = 0; k < 7; ++k) {
	  if(i==2 && k!=3) continue;
	  c.cd(26*k + i*6 + j + 1); hFieldRespInFreq2D[i][j][k]->Draw();
	  c.cd(26*k + i*6 + j + 14); hFieldRespInPhase2D[i][j][k]->Draw();
	}
      }
    }
  }
  c.SaveAs(Form("%s.pdf", canvas_name.c_str()));
}
