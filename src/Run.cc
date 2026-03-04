#include "Types.h"
#include "Utils.h"
#include "Study.h"

#include "Studies_BasicSpectra.h"
#include "Studies_QperHit.h"
#include "Studies_QvsN_2D.h"
#include "Studies_VtxXZ.h"
#include "Studies_Cone42.h"

#include "WCSimRootEvent.hh"

#include "TChain.h"
#include "TFile.h"
#include "TStyle.h"

#include <iostream>
#include <memory>
#include <vector>

int RunAnalysis(const char* fname, ConeVertexMode coneMode)
{
  gStyle->SetOptStat(0);

  TChain t("wcsimT");
  t.Add(fname);
  if (t.GetEntries() == 0) { std::cerr << "ERROR: no entries\n"; return 1; }

  std::string firstFileName = t.GetFile()->GetName();
  TFile* f = TFile::Open(firstFileName.c_str());
  if (!f || !f->IsOpen()) { std::cerr << "ERROR: could not open " << firstFileName << "\n"; return 2; }

  WCSimRootEvent* superev = new WCSimRootEvent();
  t.SetBranchAddress("wcsimrootevent", &superev);

  t.GetEntry(0);
  auto* ev0 = superev->GetTrigger(0);

  TVector3 vtx0, beamDir;
  if (!GetPrimaryVtxAndDir(ev0, vtx0, beamDir)) { std::cerr << "ERROR: no primary track id==1\n"; return 3; }

  GeoCache geo;
  if (!BuildGeoCache(f, vtx0, geo)) { std::cerr << "ERROR: BuildGeoCache failed\n"; return 4; }

  // Studies
  std::vector<std::unique_ptr<IStudy>> studies;
  studies.emplace_back(std::make_unique<BasicSpectraStudy>());
  studies.emplace_back(std::make_unique<QperHitStudy>());
  studies.emplace_back(std::make_unique<QvsN2DStudy>());
  studies.emplace_back(std::make_unique<VtxXZStudy>());
  studies.emplace_back(std::make_unique<Cone42Study>(geo, beamDir.Unit(), coneMode));

  // For summed display and per-event samples (optional)
  std::vector<double> pmt_charge_sum(geo.nPMT, 0.0);
  int savedGamma=0, savedPion=0;
  const int nSaveEach = 5;

  const long nEntries = t.GetEntries();
  long printEvery = std::max<long>(1, nEntries/100);

  long nHadronic=0, nPion=0;

  for (long iev=0; iev<nEntries; iev++) {
    if (iev % printEvery == 0) std::cout << "Running " << iev << " / " << nEntries << "\n";

    t.GetEntry(iev);
    auto* ev = superev->GetTrigger(0);

    double triggerShift=0, triggerTime=0;
    auto trigInfo = ev->GetTriggerInfo();
    if (ev->GetTriggerType()!=kTriggerNoTrig && trigInfo.size()>=3) { triggerShift=trigInfo[1]; triggerTime=trigInfo[2]; }

    EventSummary s;
    s.iev = iev;
    s.truth = ClassifyTruth(ev);
    if (s.truth.hasHadronic) nHadronic++;
    if (s.truth.hasPions) nPion++;

    double vx=0, vz=0;
    s.hasVtxXZ = GetTruthVertexXZ(ev, s.truth, vx, vz);
    if (s.hasVtxXZ) { s.vx=vx; s.vz=vz; }

    s.reco = ComputeReco(ev, geo, &pmt_charge_sum, /*hDT*/nullptr, triggerShift, triggerTime, 1e9);
    s.qPerHit = (s.reco.nDigiHits>0) ? (s.reco.totQ / s.reco.nDigiHits) : 0.0;

    for (auto& st : studies) {
        st->Fill(s);
        st->FillEvent(ev, s, triggerShift, triggerTime); // cone study uses this, others ignore it
    }

    // Sample event displays (5 no-pion and 5 pion)
    if ((savedGamma < nSaveEach && !s.truth.hasPions) || (savedPion < nSaveEach && s.truth.hasPions)) {
      std::vector<double> pmt_evt(geo.nPMT, 0.0);
      (void)ComputeReco(ev, geo, &pmt_evt, nullptr, triggerShift, triggerTime, 1e9);

      if (!s.truth.hasPions && savedGamma < nSaveEach) {
        SaveEventDisplayPDF(geo, pmt_evt, Form("event_gamma_%05d_iev%06ld.pdf", savedGamma, iev));
        savedGamma++;
      }
      if (s.truth.hasPions && savedPion < nSaveEach) {
        SaveEventDisplayPDF(geo, pmt_evt, Form("event_pion_%05d_iev%06ld.pdf", savedPion, iev));
        savedPion++;
      }
    }
  }

  // Summed display
  SaveEventDisplayPDF(geo, pmt_charge_sum, "event_display_sum.pdf");

  for (auto& st : studies) st->WritePlots();

  std::cout << "Hadronic events: " << nHadronic << " / " << nEntries << "\n";
  std::cout << "Pion events:     " << nPion << "\n";

  f->Close();
  delete f;
  delete superev;
  return 0;
}