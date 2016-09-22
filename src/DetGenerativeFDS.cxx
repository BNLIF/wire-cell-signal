#include "WireCellSignal/DetGenerativeFDS.h"
#include "TMath.h"

using namespace WireCellSignal;
using namespace WireCell;
using namespace std;

namespace {
  double integral(double sigma, double z0, double z1, double z2){
    double result=0;
    result += 0.5 * TMath::Erf((z2-z0)/sqrt(2.)/sigma) 
      - 0.5 *TMath::Erf((z1-z0)/sqrt(2.)/sigma);
    return result;
  }
}

WireCellSignal::DetGenerativeFDS::DetGenerativeFDS(const Depositor& dep, const DetectorGDS& gds, 
				   int bins_per_frame1, int nframes_total,
				   float bin_drift_distance, float unit_dis,
				   float longest_drift)
    : dep(dep)
    , det_gds(gds)
    , max_frames(nframes_total)
    , bin_drift_distance(bin_drift_distance)
    , unit_dis(unit_dis)
    , max_drift(longest_drift) // might have to consider the possibility that the maximum drift distance is not the same for both sides
{
  bins_per_frame = bins_per_frame1;
  tot_charge = 0;
  collected_charge = 0;
  const int napas = gds.napa(0);
  sideband_charge[0] = new double[napas];
  sideband_charge[1] = new double[napas];
  for (int i = 0; i < napas; ++i) {
    sideband_charge[0][i] = 0;
    sideband_charge[1][i] = 0;
  }
  for (int cryo = 0; cryo<det_gds.ncryos(); ++cryo) {
    for (int apa = 0; apa<det_gds.napa(cryo); ++apa) {
      const WrappedGDS* apa_gds = det_gds.get_apaGDS(cryo, apa);	
      for (int i = 0; i < 3; ++i) {
	WirePlaneType_t plane = static_cast<WirePlaneType_t>(i);
	for (int face = 0; face<2; ++face) {
	  int nwires = apa_gds->wires_in_plane(face,plane).size();
	  for (int j = 0; j < 10; ++j) {	    
	    hraw[cryo][apa][i][face][j] = new TH2F(Form("hraw_%d_%d_%d_%d_%d", cryo, apa, i, face, j),
						   Form("hraw_%d_%d_%d_%d_%d", cryo, apa, i, face, j),
						   nwires, 0 , nwires,
						   bins_per_frame, 0, bins_per_frame);
	    /*
	      hconv[cryo][apa][i][face][j] = new TH2D(Form("hconv_%d_%d_%d_%d_%d", cryo, apa, i, face, j),
						   Form("hconv_%d_%d_%d_%d_%d", cryo, apa, i, face, j),
						   nwires, 0 , nwires,
						   bins_per_frame, 0, bins_per_frame);
	    */
	  }
	}
      }
    }
  }
}

WireCellSignal::DetGenerativeFDS::~DetGenerativeFDS()
{
}

int WireCellSignal::DetGenerativeFDS::size() const { return max_frames; }

WireCell::SimTruthSelection WireCellSignal::DetGenerativeFDS::truth() const
{
    WireCell::SimTruthSelection ret;

    for (auto it = simtruth.begin(); it != simtruth.end(); ++it) {
	ret.push_back(& (*it));
    }
    return ret;
}

int WireCellSignal::DetGenerativeFDS::jump(int frame_number)
{  
    if (frame_number < 0) {
	return -1;
    }
    if (max_frames >=0 && frame_number >= max_frames) {
	return -1;
    }

    if (frame.index == frame_number) {
	return frame_number;
    }

    frame.clear();
    simtruth.clear();

    const PointValueVector& hits = dep.depositions(frame_number);
    const std::vector<int>& timeoffsets = dep.timeoffset();
    //const double eCharge = 1.602177e-4; // charge of an electron in fC
    
    size_t nhits = hits.size();
    //std::cout<<"nhits = "<<nhits<<std::endl;
    if (!nhits) {
	frame.index = frame_number; 
	return frame.index;
    }

    //std::cout<<"frame.index != frame_number"<<std::endl;
    
    typedef map<int, int> TraceIndexMap; // channel->index into traces;
    TraceIndexMap tim;		// keep tabs on what channels we've seen already

    std::pair<double, double> xmm;
    
    //std::cout << nhits << std::endl;

    const WrappedGDS *apa_gds = 0;

    for (size_t ind=0; ind<nhits; ++ind) {

      //if (ind<1000) std::cout<<ind<<std::endl;
      //std::cout<<ind<<std::endl;
      //if (ind > 1000) break;

      const Point& pt = hits[ind].first;
      double charge = hits[ind].second;//*eCharge;
      tot_charge += charge;      
      if (ind%1000==0) std::cout<<ind<<", charge = "<<charge<<std::endl;
           
      // decide which apa the charge deposit belongs to
      short which_cryo = det_gds.in_which_cryo(pt);
      short which_apa = det_gds.in_which_apa(pt);
      apa_gds = det_gds.get_apaGDS(det_gds.in_which_cryo(pt), det_gds.in_which_apa(pt));
      if (apa_gds==NULL) {
	  // out side the detector boundary parallel to drift direction 
	std::cout<<"exit because apa_gds = NULL"<<std::endl; 
	continue;
      }
      xmm = apa_gds->minmax(0); 
      
      //	std::cout << pt.x/units::cm << " " << xmm.first/units::cm << " " << xmm.second/units::cm << std::endl;
      if (pt.x>xmm.first && pt.x<xmm.second) continue;// space in between wires in an APA are treated as dead region
      
      int face = 0; // "A face -x"
      float drift_dist;
      
      if (TMath::Abs(pt.x-xmm.first) > TMath::Abs(pt.x-xmm.second)) {
	drift_dist = TMath::Abs(pt.x-xmm.second);
	face = 1; //"B face +x"
      }else{
	drift_dist = TMath::Abs(pt.x-xmm.first);
      }
      if (drift_dist > max_drift) continue; // out of TPC boundary along drift direction

      // float drift_dist = TMath::Abs(pt.x-xmm.first) < TMath::Abs(pt.x-xmm.second) ? TMath::Abs(pt.x-xmm.first) : TMath::Abs(pt.x-xmm.second);
      int tbin = int(drift_dist/bin_drift_distance);
      tbin = TMath::Abs(tbin);
      int offset;
      //std::cout << timeoffsets.size() << std::endl;
      if (timeoffsets.size()==0){
	offset = 0;
      }else{
	offset = timeoffsets[ind];
      }
      //std::cout << tbin << " " << pt.x/units::cm << " " << xmm.second/units::cm << " " << offset << std::endl;
      
      // adding in the diffusion
      // assuming the velocity is 1.6 mm/us
      float drift_time = drift_dist/(unit_dis*units::millimeter); // us      
      
      // can not handle negative drift time yet .... 
      if (drift_time < 0) {
	//std::cout<<"exit because drift time < 0"<<std::endl; 
	continue;
      }
      
      float DL = 5.3; //cm^2/s
      float DT = 12.8; //cm^2/s
      //float DL = 0.1;
      //float DT = 0.1;
      float sigmaL = sqrt(2.*DL*drift_time*1e-6) * units::cm;
      float sigmaT = sqrt(2.*DT*drift_time*1e-6) * units::cm;

      // do longitudinal array first
      int ntbin = sigmaL*3/bin_drift_distance + 1; // +- ntinb 3sigma
      std::vector<int> long_tbin;
      std::vector<double> long_integral;
      //push middle bin first
      long_tbin.push_back(tbin);
      long_integral.push_back(integral(sigmaL,drift_dist,tbin*bin_drift_distance,(tbin+1)*bin_drift_distance));
            
      for (int kk =0; kk!=ntbin;kk++){
	int tt = tbin - kk - 1;
	if (tt >0){
	  long_tbin.push_back(tt);
	  long_integral.push_back(integral(sigmaL,drift_dist,tt*bin_drift_distance,(tt+1)*bin_drift_distance));
	}
	
	tt = tbin + kk + 1;
	if (tt < bins_per_frame){
	  long_tbin.push_back(tt);
	  long_integral.push_back(integral(sigmaL,drift_dist,tt*bin_drift_distance,(tt+1)*bin_drift_distance));
	}
      }
            
      if (tbin >= bins_per_frame) {
	//cerr << "GenerativeFDS: drop: tbin=" << tbin << " >= " << bins_per_frame << endl;
	continue;
      }
      
      WireCell::SimTruth st(pt.x, pt.y, pt.z, charge, tbin, simtruth.size());
      simtruth.insert(st);

      for (int iplane=0; iplane < 3; ++iplane) {
	WirePlaneType_t plane = static_cast<WirePlaneType_t>(iplane); // annoying
	  // only look up wires from the correct face
	const GeomWire* wire = apa_gds->closest(pt, plane, face);
	const GeomWireSelection& wires_in_plane = apa_gds->wires_in_plane(face, plane);
	const GeomWire *wirefirst=wires_in_plane.at(0);
	const GeomWire *wirelast=wires_in_plane.at(wires_in_plane.size()-1);
	//const GeomWire *wire50=wires_in_plane.at(50);
	double wirestart = apa_gds->wire_dist(*wirefirst) - apa_gds->pitch(0,plane)/2.;
	double wireend  = apa_gds->wire_dist(*wirelast) + apa_gds->pitch(0,plane)/2.;
	//double randm = apa_gds->wire_dist(*wire50);
	//std::cout<<"plane = "<<iplane<<", size = "<<wires_in_plane.size()<<", zerox = "<<zerox<<", end = "<<maxx<<", random = "<<randm<<std::endl;
	if (wire!=0){
	  int chid = wire->channel();
	  int windex = wire->index();
	  //int wirebin = (int)(find(wires_in_plane.begin(), wires_in_plane.end(), wire) - wires_in_plane.begin());
	  //std::cout<<"wirebin = "<<wirebin<<std::endl;//////////////////
	  // start to do the transverse diffusion here ...
	  double pitch = apa_gds->pitch(0,plane);
	  double angle = apa_gds->angle(plane);
	  double increment = pitch/10.;
	  double dist = apa_gds->wire_dist(pt,plane,wire->face());
	  double wiredist = apa_gds->wire_dist(*wire);
	  int nwbin = sigmaT*3 / increment;

	  std::vector<double> trans_integral;
	  std::vector<double> trans_dist;
	  for (int kk = -nwbin; kk <= nwbin; ++kk) {
	    trans_integral.push_back( integral(sigmaT, dist, wiredist-(double)kk*increment-increment/2., wiredist-(double)kk*increment+increment/2.) );
	    trans_dist.push_back(wiredist-(double)kk*increment);
	  }

	  std::vector<double> trans_charge;
	  std::vector<int> trans_index1; // point index at which field response is simulated
	  std::vector<int> trans_wire_bin;
	  for (std::size_t kk = 0; kk != trans_dist.size(); ++kk) {
	    double dis = trans_dist.at(kk);
	    if (dis<wirestart || dis>wireend) continue;
	    const GeomWire *wire1 = apa_gds->closest(dis, plane, face);
	    if (wire1) {
	      int wire1bin = trunc((dis-wirestart)/pitch);
	      int index1 = trunc((dis-wirestart-wire1bin*pitch)/pitch*10);
	      if (index1>9 || index1<0) {
		std::cerr << "FIXME!! index1 = "<<index1<<" when it shouldn't be <0 or >10!" << std::endl;
		continue;
	      }
	      trans_wire_bin.push_back(wire1bin);
	      trans_index1.push_back(index1);
	      trans_charge.push_back(trans_integral.at(kk));
	    }
	  }
	  /*
	  GeomWireSelection trans_wires;
	  std::vector<double> trans_charge;
	  std::vector<int> trans_index1; // point index at which field response is simulated
	  std::vector<int> trans_index2; // wire index associated with a field response
	  std::vector<int> trans_wire_bin;
	  for (std::size_t kk = 0; kk != trans_dist.size(); ++kk) {
	    for (int ii = -3; ii < 4; ++ii) {
	      double dis = trans_dist.at(kk) + (double)ii * pitch;
	      const GeomWire *wire1 = apa_gds->closest(dis, plane, face);
	      if (wire1) {
		int wire1bin = (int)(find(wires_in_plane.begin(), wires_in_plane.end(), wire1) - wires_in_plane.begin());
		double wire1_dis = apa_gds->wire_dist(*wire1);		
		int index1 = round((dis-wire1_dis)/increment) + 4;
		if (index1>9 || index1<0) {
		  //std::cerr << "FIXME!! index1 = "<<index1<<" when it shouldn't be <0 or >10!" << std::endl;
		  continue;
		}
		trans_wire_bin.push_back(wire1bin);
		trans_index1.push_back(index1);
		if (TMath::Abs(trans_dist.at(kk) - wire1_dis) > TMath::Abs(pitch * ii))
		  trans_index2.push_back(ii);
		else trans_index2.push_back(-1*ii);
		trans_wires.push_back(wire1);
		trans_charge.push_back(trans_integral.at(kk));
	      }
	    }
	  }
	  */	  
	  //GeomWireSelection allwires;
	  std::vector<int> alltime;
	  std::vector<float> allcharge;
	  std::vector<int> allindex1;
	  //std::vector<int> allindex2;
	  std::vector<int> allwirebin;
	  float sum_charge = 0;
	  for (int qt = 0; qt!= long_tbin.size(); qt++){
	    for (int qw = 0; qw!= trans_charge.size(); qw++){
	      alltime.push_back(long_tbin.at(qt) + offset);
	      //allwires.push_back(trans_wires.at(qw));
	      float tcharge = long_integral.at(qt) * 
		trans_charge.at(qw);
	      allcharge.push_back(tcharge);
	      allwirebin.push_back(trans_wire_bin.at(qw));
	      allindex1.push_back(trans_index1.at(qw));	      
	      //allindex2.push_back(trans_index2.at(qw));
	      //if (trans_index2.at(qw)==0) sum_charge += tcharge; // amounts to collected charge
	      sum_charge += tcharge;
	    }
	  }
	  // do normalization and fill TH2D
	  for (int qx = 0; qx!=allcharge.size(); ++qx) {
	    int tbin3 = alltime.at(qx);
	    if (tbin3 >=0 && tbin3 < bins_per_frame){
	      float charge3 = 0;
	      if (sum_charge>0) charge3 = allcharge.at(qx)/sum_charge*charge;
	      //const GeomWire* wire3 = allwires.at(qx);
	      int wire_bin = allwirebin.at(qx);
	      //std::cout<<wire_bin<<", "<<charge3<<", "<<alltime.at(qx)<<std::endl;
	      hraw[which_cryo][which_apa][iplane][face][allindex1.at(qx)]->Fill(wire_bin,alltime.at(qx), charge3);
	    }
	  }
	  /*	  
	  //do normalization ... 
	  for (int qx = 0; qx!=allcharge.size();qx++){
	    const GeomWire* wire3 = allwires.at(qx);
	    int chid3 = wire3->channel();
	    int tbin3 = alltime.at(qx);	      
	    if (tbin3 >=0 && tbin3 < bins_per_frame){		
	      float charge3 = 0;
	      if (sum_charge>0) charge3 = allcharge.at(qx)/sum_charge*charge;
	      //std::cout<<"total charge = "<<charge<<", charge3 = "<<charge3<<std::endl;
	      TraceIndexMap::iterator it = tim.find(chid3);
	      int trace_index = frame.traces.size(); // if new
	      if (it == tim.end()) {
		Trace t;
		t.chid = chid3;
		t.tbin = 0; // Need change???
		t.hCharge = new TH1F("sig","",bins_per_frame, 0, bins_per_frame); // make sure bins_per_frame is consistent with the definition in electronics configuration
		t.channel_length = apa_gds->get_channel_length(chid3);
		tim[chid3] = frame.traces.size();
		frame.traces.push_back(t);
	      }else {		// already seen
		trace_index = it->second;
	      }
	      Trace& trace = frame.traces[trace_index];	      
	      // finally			       
	      //trace.charge[tbin3] += charge3;
	      hraw[plane][allindex2.at(qx)+5]->Fill(chid3, tbin3, charge3);
	      //trace.hCharge->Add(fRsp->GetResponseFunction(plane, allindex1.at(qx), allindex2.at(qx), tbin3, charge3));	      
	      if (plane==2) {
		collected_charge += charge3;
		if (drift_dist < 5.*units::cm) sideband_charge[face][which_apa]+= charge3;
		//std::cout<<"charge3 = "<<charge3<<", collected_charge = "<<collected_charge<<", sum_charge = "<<sum_charge<<", charge = "<<charge<<std::endl;
	      }
	    }
	  }
	  */
	}
	
      }
    }
    /*
    TFile *tmpfile = new TFile("real_signal.root","recreate");
    return 0;

    for (int i=0; i<det_gds.ncryos(); ++i) {
      for (int j=0; j<det_gds.napa(i); ++j) {
	for (int k=0; k<3; ++k) {
	  for (int l=0; l<2; ++l) {
	    for (int m=0; m<10; ++m) {
	      if (hraw[i][j][k][l][0]->GetEntries()>0) {
		hraw[i][j][k][l][0]->Add(hraw[i][j][k][l][m]);
	      }
	      if (k==0) {
		for (int bx=0; bx<hraw[i][j][k][l][0]->GetNbinsX(); ++bx) {
		  for (int by=0; by<hraw[i][j][k][l][0]->GetNbinsY(); ++by) {
		    if (bx<450) hraw[i][j][k][l][m]->SetBinContent(bx+1, by+1,hraw[i][j][k][l][m]->GetBinContent(bx+1, by+1)*10);
		  }
		}
	      }
	      if (k==2) {
		for (int bx=0; bx<hraw[i][j][k][l][0]->GetNbinsX(); ++bx) {
		  for (int by=0; by<hraw[i][j][k][l][0]->GetNbinsY(); ++by) {
		    if (bx<55 && by<3440) hraw[i][j][k][l][m]->SetBinContent(bx+1, by+1, hraw[i][j][k][l][m]->GetBinContent(bx+1, by+1)*10);
		  }
		}
	      }

	    }
	    //tmpfile->cd();
	    //hraw[i][j][k][l][0]->Write();
	  }
	}
      }
    }   
    */
    
    // start 2D convolution
    std::cout<<"\nprepare 2D convolution\n"<<std::endl;
    TStopwatch *sw = new TStopwatch();
    const int dim = bins_per_frame;
    // do 2D convolution
    std::cout<<"start actual convolution\n"<<std::endl;    
    TH1F *htemp = new TH1F("htemp","htemp",bins_per_frame,0,bins_per_frame);

    for (int i=0; i<det_gds.ncryos(); ++i) {
      for (int j=0; j<det_gds.napa(i); ++j) {
	for (int k=0; k<3; ++k) {
	  for (int m=0; m<10; ++m) {
	    sw->Reset();
	    sw->Start();
	    double **rho_res, **phi_res;
	    rho_res = new double*[7];
	    phi_res = new double*[7];
	    for (int w=0; w<7; ++w) {
	      rho_res[w] = new double[dim];
	      phi_res[w] = new double[dim];
	    }
	    //double rho_res[7][dim], phi_res[7][dim];
	    for (int w=0; w!=7; ++w) {
	      int wireid = w-3;
	      if (m-4<0) wireid = 3-w;
	      TH1 *tmp = fRsp->GetResponseFunction(k, abs(m-4), wireid);
	      TVirtualFFT::SetTransform(0);
	      TH1 *hph = tmp->FFT(0,"PH");
	      TH1 *hmag = tmp->FFT(0,"MAG");
	      /*
	      if (k==2 && w!=3) {
		hph->Reset();
		hmag->Reset();
	      }
	      */
	      for (int t = 0; t != bins_per_frame; ++t) {
		rho_res[w][t] = hmag->GetBinContent(t+1);
		phi_res[w][t] = hph->GetBinContent(t+1);
	      }
	      delete hph;
	      delete hmag;
	    }
	    
	    for (int l=0; l<2; ++l) {	      
	      if (hraw[i][j][k][l][m]->GetEntries()>0) {
		
		const int xdim = hraw[i][j][k][l][m]->GetNbinsX(); // wires
		std::cout<<"# of wires is "<<xdim<<std::endl;
		double **rho_v, **phi_v;
		double **result_re, **result_im;		
		rho_v = new double*[xdim];
		phi_v = new double*[xdim];
		result_re = new double*[xdim];
		result_im = new double*[xdim];
		for (int w=0; w!=xdim; ++w) {
		  rho_v[w] = new double[bins_per_frame];
		  phi_v[w] = new double[bins_per_frame];
		  result_re[w] = new double[bins_per_frame];
		  result_im[w] = new double[bins_per_frame];
		}
		for (int w=0; w!=xdim; ++w){
		  htemp->Reset();
		  for (int t=0; t!=bins_per_frame; ++t) {
		    htemp->SetBinContent(t, hraw[i][j][k][l][m]->GetBinContent(w+1, t+1));
		  }
		  TH1 *hm = htemp->FFT(0,"MAG");
		  TH1 *hp = htemp->FFT(0,"PH");
		  for (Int_t t=0;t!=bins_per_frame;t++){
		    rho_v[w][t] = hm->GetBinContent(t+1);
		    phi_v[w][t] = hp->GetBinContent(t+1);
		  }
		  delete hm;
		  delete hp;
		}
		std::cout<<"fft-ed signal"<<std::endl;

		//double result_re[xdim][dim], result_im[xdim][dim];
		double *value_re; value_re = new double[dim];
		double *value_im; value_im = new double[dim];
		double *resp_re;  resp_re = new double[xdim];
		double *resp_im;  resp_im = new double[xdim];
		//double value_re[dim], value_im[dim];
		//double resp_re[xdim], resp_im[xdim];
		for (int t=0; t!=bins_per_frame; ++t) {
		  for (int w=0; w!=xdim; ++w) {
		    value_re[w] = rho_v[w][t]*cos(phi_v[w][t]);
		    value_im[w] = rho_v[w][t]*sin(phi_v[w][t]);
		    if (w<7) {
		      resp_re[w] = rho_res[w][t]*cos(phi_res[w][t]);
		      resp_im[w] = rho_res[w][t]*sin(phi_res[w][t]);
		      //std::cout<<value_re[w]<<" "<<value_im[w]<<" "<<resp_re[w]<<" "<<resp_im[w]<<std::endl;
		    } else {
		      resp_re[w] = 0;
		      resp_im[w] = 0;
		    }
		  }
		  
		  int m = xdim;
		  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&m,"C2CFORWARD M K");
		  ifft->SetPointsComplex(value_re,value_im);
		  ifft->Transform();
		  Double_t temps_re[xdim],temps_im[xdim];
		  ifft->GetPointsComplex(temps_re,temps_im);

		  ifft->SetPointsComplex(resp_re,resp_im);
		  ifft->Transform();
		  Double_t tempw_re[xdim],tempw_im[xdim];
		  ifft->GetPointsComplex(tempw_re,tempw_im);
		  delete ifft;
		  
		  Double_t tempr_re[xdim],tempr_im[xdim];
		  for (int w=0; w!=xdim; ++w) {
		    //std::cout<<tempw_re[w]<<" "<<tempw_im[w]<<" "<<temps_re[w]<<" "<<temps_im[w]<<std::endl;
		    tempr_re[w] = (temps_re[w]*tempw_re[w] - temps_im[w]*tempw_im[w])/xdim;
		    tempr_im[w] = (temps_im[w]*tempw_re[w] + temps_re[w]*tempw_im[w])/xdim;
		    //std::cout<<tempr_re[w]<<" "<<tempr_im[w]<<std::endl;
		  }
		  
		  TVirtualFFT *ifft3 = TVirtualFFT::FFT(1,&m,"C2CBACKWARD M K");
		  ifft3->SetPointsComplex(tempr_re,tempr_im);
		  ifft3->Transform();
		  Double_t temp3_re[xdim],temp3_im[xdim];
		  ifft3->GetPointsComplex(temp3_re,temp3_im);
		  delete ifft3;
		  for (Int_t w=0; w!=xdim; ++w){
		    Int_t shift = w - 3;
		    if (shift <0) shift += xdim;
		    result_re[w][t] = temp3_re[shift]/bins_per_frame;
		    result_im[w][t] = temp3_im[shift]/bins_per_frame;
		    //if (result_re[w][t]>1e-5 || result_im[w][t]>1e-5) std::cout<<result_re[w][t]<<" "<<result_im[w][t]<<std::endl;
		  }
		  
		}
		delete[] value_re;
		delete[] value_im;
		delete[] resp_re;
		delete[] resp_im;
		for (int w=0; w!=xdim; ++w) {
		  delete[] rho_v[w];
		  delete[] phi_v[w];
		}
		delete[] rho_v;
		delete[] phi_v;
		std::cout<<"convoluted in frequency space"<<std::endl;
		
		double temp_re[dim],temp_im[dim];
		std::cout<<"start fft back to time space"<<std::endl;
		hraw[i][j][k][l][m]->Reset();
		std::cout<<"reset hraw[i][j][k][l][m]"<<std::endl;
		for (int w=0; w!=xdim; ++w) {
		  int n = bins_per_frame;
		  TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
		  for (int t=0; t!=bins_per_frame; ++t){
		    temp_re[t] = result_re[w][t];
		    temp_im[t] = result_im[w][t];
		    //std::cout<<temp_re[t]<<" "<<temp_im[t]<<std::endl;
		  }
		  ifft2->SetPointsComplex(temp_re,temp_im);
		  ifft2->Transform();		  
		  TH1 *fb = 0;
		  fb = TH1::TransformHisto(ifft2,fb,"Re");
		  //std::cout<<"got fb "<<std::endl;//fb->GetBinContent(2)<<std::endl;
		  for (int t=0; t!=bins_per_frame; ++t){
		    hraw[i][j][k][l][m]->SetBinContent(w+1, t+1, fb->GetBinContent(t+1));
		    //std::cout<<t<<std::endl;
		  }
		  delete ifft2;
		  delete fb;
		  //std::cout<<"end"<<std::endl;
		}
		for (int w=0; w!=xdim; w++) {
		  delete[] result_re[w];
		  delete[] result_im[w];
		}
		delete[] result_re;
		delete[] result_im;
		std::cout<<"end of hraw[i][j][k][l][m]->GetEntries()>0\n"<<std::endl;
	      } // end of if hraw[i][j][k][l][m]->GetEntries()>0
	    }
	    for (int w=0; w!=7; ++w) {
	      delete[] rho_res[w];
	      delete[] phi_res[w];
	    }
	    delete[] rho_res;
	    delete[] phi_res;
	    
	    sw->Stop();
	    sw->Print("u");
	  }
	}
      }
    }

    for (int i=0; i<det_gds.ncryos(); ++i) {
      for (int j=0; j<det_gds.napa(i); ++j) {
	for (int k=0; k<3; ++k) {
	  for (int l=0; l<2; ++l) {
	    for (int m=1; m<10; ++m) {
	      if (hraw[i][j][k][l][m]->GetEntries()) {
		hraw[i][j][k][l][0]->Add(hraw[i][j][k][l][m]);
		delete hraw[i][j][k][l][m];
	      }
	    }
	  }
	}
      }
    }
    
    std::cout<<"before saving"<<std::endl;
    TFile f("tmp.root","recreate");
    f.cd();
    for (int i=0; i<det_gds.ncryos(); ++i) {
      for (int j=0; j<det_gds.napa(i); ++j) {
	const WrappedGDS *apa_gds = det_gds.get_apaGDS(i,j);
	for (int k=0; k<3; ++k) {
	  for (int l=0; l<2; ++l) {
	    if (hraw[i][j][k][l][0]!=NULL) hraw[i][j][k][l][0]->Write();
	    const GeomWireSelection& wires_in_plane = apa_gds->wires_in_plane(l, (WirePlaneType_t)k);
	    for (int bx=0; bx<hraw[i][j][k][l][0]->GetNbinsX(); ++bx) {
	      const GeomWire* wire = wires_in_plane.at(bx);
	      int chid = wire->channel();
	      TraceIndexMap::iterator it = tim.find(chid);
	      int trace_index = frame.traces.size();
	      if (it == tim.end()) {
		Trace t;
		t.chid = chid;
		t.tbin = 0;
		t.hCharge = new TH1F(Form("rs%d",chid), "", bins_per_frame, 0, bins_per_frame);
		t.channel_length = apa_gds->get_channel_length(chid);
		tim[chid] = frame.traces.size();
		frame.traces.push_back(t);
	      } else {
		trace_index = it->second;
	      }
	      Trace& trace = frame.traces[trace_index];
	      for (int by=0; by<hraw[i][j][k][l][0]->GetNbinsY(); ++by) {
		trace.hCharge->SetBinContent(by+1, trace.hCharge->GetBinContent(by+1) + hraw[i][j][k][l][0]->GetBinContent(bx+1,by+1));
	      }
	    }
	  }
	}
      }
    }
    f.Close();
    
    frame.index = frame_number; 
    return frame.index;
}

