#pragma once
#include "Study.h"
#include "Types.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"

class WCSimRootTrigger;

enum class ConeVertexMode { Nominal, Truth };

struct ConeStudy : IStudy {
  const GeoCache& geo;
  TVector3 beam_hat;          // unit vector (beam axis)
  ConeVertexMode vtxMode;
  TVector3 vtx_nom;           // cm

  double cone_deg;            // half-angle in degrees
  double cosCone;             // cos(cone_deg)

  // histograms (fractions)
  TH1D *h_fNout_pion, *h_fNout_nopion;
  TH1D *h_fQout_pion, *h_fQout_nopion;
  TH1D *h_fNout_emshower_pion; // subset of pion events without hadronic activity
  TH1D *h_fQout_emshower_pion; // subset of pion events without hadronic activity

  // same-side
  TH2D *h2_Qin_Nin_pion,   *h2_Qin_Nin_nopion;
  TH2D *h2_Qout_Nout_pion, *h2_Qout_Nout_nopion;

  // cross
  TH2D *h2_Qin_Nout_pion,  *h2_Qin_Nout_nopion;
  TH2D *h2_Qout_Nin_pion,  *h2_Qout_Nin_nopion;

  // fractional observables
  TH2D *h2_fQ_fN_pion;
  TH2D *h2_fQ_fN_nopion;

  ConeStudy(const GeoCache& g,
            const TVector3& beam_unit,
            ConeVertexMode mode,
            double coneDeg);

  void Fill(const EventSummary& /*e*/) override {}
  void FillEvent(WCSimRootTrigger* ev,
                 const EventSummary& e,
                 double triggerShift,
                 double triggerTime) override;

  void WritePlots() override;

private:
  TVector3 PickVertex(WCSimRootTrigger* ev, const EventSummary& e) const;
  std::string Tag() const; // e.g. "Cone42" or "Cone60"
};
