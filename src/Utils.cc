#include "Utils.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <cmath>

bool GetPrimaryVtxAndDir(WCSimRootTrigger* ev, TVector3& vtx, TVector3& dir) {
  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (trk && trk->GetId() == 1) {
      vtx = TVector3(trk->GetStart(0), trk->GetStart(1), trk->GetStart(2));
      dir = TVector3(trk->GetDir(0),   trk->GetDir(1),   trk->GetDir(2));
      return true;
    }
  }
  return false;
}

bool BuildGeoCache(TFile* f, const TVector3& vtx, GeoCache& geo) {
  auto* geotree = (TTree*) f->Get("wcsimGeoT");
  if (!geotree) { std::cerr << "ERROR: could not find wcsimGeoT\n"; return false; }

  WCSimRootGeom* g = nullptr;
  geotree->SetBranchAddress("wcsimrootgeom", &g);
  if (geotree->GetEntries() <= 0) { std::cerr << "ERROR: wcsimGeoT has 0 entries\n"; return false; }
  geotree->GetEntry(0);

  geo.nPMT = g->GetWCNumPMT();
  geo.vg_cm_per_ns = 2.20027795333758801e8 * 100 / 1.e9;

  geo.pmt_pos.resize(geo.nPMT);
  geo.pmt_tof.resize(geo.nPMT);
  geo.pmt_cylLoc.resize(geo.nPMT);
  geo.mpmtNo.resize(geo.nPMT);
  geo.mpmt_pmtNo.resize(geo.nPMT);
  geo.eventXY.resize(geo.nPMT);

  geo.max_r = 0.0;
  geo.max_z = 0.0;

  for (int i=0; i<geo.nPMT; i++) {
    WCSimRootPMT pmt = g->GetPMT(i);

    double px = pmt.GetPosition(0);
    double py = pmt.GetPosition(1);
    double pz = pmt.GetPosition(2);

    geo.mpmtNo[i]      = pmt.GetmPMTNo();
    geo.mpmt_pmtNo[i]  = pmt.GetmPMT_PMTNo();
    geo.pmt_cylLoc[i]  = pmt.GetCylLoc();

    TVector3 pmtpos(px, py, pz);
    geo.pmt_pos[i] = pmtpos;
    geo.pmt_tof[i] = (pmtpos - vtx).Mag() / geo.vg_cm_per_ns;

    geo.max_z = std::max(geo.max_z, std::fabs(py));

    double r = std::sqrt(px*px + pz*pz);
    if (std::fabs(pmt.GetOrientation(1)) > 0.5) geo.max_r = std::max(geo.max_r, r);
  }

  geo.barrelCut = geo.max_z - 10.0;

  for (int i=0; i<geo.nPMT; i++) {
    double y = -geo.pmt_pos[i].x();
    double x =  geo.pmt_pos[i].z();
    double z =  geo.pmt_pos[i].y();

    std::vector<double> xy(2,0.0);
    if (std::fabs(z) < geo.barrelCut) {
      double th = std::atan2(y,x);
      xy[0] = -geo.max_r * th;
      xy[1] = z;
    } else if (z > geo.barrelCut) {
      xy[0] = -y;
      xy[1] = geo.max_z + geo.max_r - x;
    } else {
      xy[0] = -y;
      xy[1] = -geo.max_z - geo.max_r + x;
    }
    geo.eventXY[i] = xy;
  }

  std::cout << "Loaded geometry: nPMT=" << geo.nPMT
            << " max_r=" << geo.max_r
            << " max_z=" << geo.max_z << "\n";
  return true;
}

EventTruth ClassifyTruth(WCSimRootTrigger* ev) {
  EventTruth t;
  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (!trk) continue;
    int pdg = trk->GetIpnu();
    int parent = trk->GetParentId();

    if (std::abs(pdg)==2212 || std::abs(pdg)==2112 || std::abs(pdg)==211 || pdg==111) t.hasHadronic = true;
    if (pdg==111 && parent==1) { t.hasPions=true; t.hasPi0=true; }
    if ((pdg==211 || pdg==-211) && parent==1) { t.hasPions=true; t.hasPiCh=true; }
  }
  return t;
}

EventReco ComputeReco(WCSimRootTrigger* ev, const GeoCache& geo,
                      std::vector<double>* pmt_charge_accum,
                      TH1D* hDigiTimeMinusTOF,
                      double triggerShift, double triggerTime,
                      double timeCutNs)
{
  EventReco r;
  const int nDigi = ev->GetNcherenkovdigihits();
  double tfirst = 1e9;
  double tavg = 0.0;
  for (int i=0; i<nDigi; i++) {
    auto* dh = (WCSimRootCherenkovDigiHit*) ev->GetCherenkovDigiHits()->At(i);
    if (!dh) continue;
    int tube = dh->GetTubeId() - 1;
    if (tube < 0 || tube >= geo.nPMT) continue;

    double q = dh->GetQ();
    double t = dh->GetT() + triggerTime - triggerShift;
    if (t > timeCutNs) continue;

    r.nDigiHits++;
    r.totQ += q;
    if (pmt_charge_accum) (*pmt_charge_accum)[tube] += q;
    if (hDigiTimeMinusTOF) hDigiTimeMinusTOF->Fill(t - geo.pmt_tof[tube], q);
    if (t < tfirst) tfirst = t;
    tavg += t;
  }
  r.t = (r.nDigiHits > 0) ? tfirst : 0.0;
  // need to be careful not to divide by zero here; tavg is only meaningful if nDigiHits > 0
  r.tavg = (r.nDigiHits > 0) ? (tavg / r.nDigiHits) : 0.0;  
  
  return r;
}

bool GetTruthVertexXZ(WCSimRootTrigger* ev, const EventTruth& truth, double& vx, double& vz) {
  if (truth.hasPions) {
    for (int i=0; i<ev->GetNtrack(); i++) {
      auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
      if (!trk) continue;
      int pdg = trk->GetIpnu();
      if (pdg==111 || pdg==211 || pdg==-211) { vx=trk->GetStart(0); vz=trk->GetStart(2); return true; }
    }
  }
  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (!trk) continue;
    int pdg = trk->GetIpnu();
    if (pdg==11 || pdg==-11) { vx=trk->GetStart(0); vz=trk->GetStart(2); return true; }
  }
  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (trk && trk->GetId()==1) { vx=trk->GetStart(0); vz=trk->GetStart(2); return true; }
  }
  return false;
}

void NormalizeToUnitArea(TH1* h) {
  double I = h->Integral(0, h->GetNbinsX()+1);
  if (I > 0) h->Scale(1.0/I);
}

void SaveEventDisplayPDF(const GeoCache& geo,
                         const std::vector<double>& pmt_charge,
                         const std::string& outname)
{
  TH2D h(Form("h_%s", outname.c_str()), "Event display;X;Y",
         250, -TMath::Pi()*geo.max_r,  TMath::Pi()*geo.max_r,
         250, -geo.max_z-2*geo.max_r, geo.max_z+2*geo.max_r);
  h.SetDirectory(nullptr);  // IMPORTANT: avoid gDirectory ownership

  for (int ipmt=0; ipmt<geo.nPMT; ipmt++) h.Fill(geo.eventXY[ipmt][0], geo.eventXY[ipmt][1], pmt_charge[ipmt]);

  TCanvas c(Form("c_%s", outname.c_str()), "c", 900, 700);
  h.Draw("colz");
  c.SaveAs(outname.c_str());
}

bool GetTruthVertexXYZ(WCSimRootTrigger* ev, const EventTruth& truth,
                       double& vx, double& vy, double& vz)
{
  if (truth.hasPions) {
    for (int i=0; i<ev->GetNtrack(); i++) {
      auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
      if (!trk) continue;
      int pdg = trk->GetIpnu();
      if (pdg==111 || pdg==211 || pdg==-211) {
        vx = trk->GetStart(0); vy = trk->GetStart(1); vz = trk->GetStart(2);
        return true;
      }
    }
  }

  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (!trk) continue;
    int pdg = trk->GetIpnu();
    if (pdg==11 || pdg==-11) {
      vx = trk->GetStart(0); vy = trk->GetStart(1); vz = trk->GetStart(2);
      return true;
    }
  }

  for (int i=0; i<ev->GetNtrack(); i++) {
    auto* trk = (WCSimRootTrack*) ev->GetTracks()->At(i);
    if (trk && trk->GetId()==1) {
      vx = trk->GetStart(0); vy = trk->GetStart(1); vz = trk->GetStart(2);
      return true;
    }
  }
  return false;
}
