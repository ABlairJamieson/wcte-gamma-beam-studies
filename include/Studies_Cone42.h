#pragma once
#include "Study.h"
#include "Types.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>

class WCSimRootTrigger; 

enum class ConeVertexMode { Nominal, Truth };

struct Cone42Study : IStudy {
  const GeoCache& geo;
  TVector3 beam_hat;              // unit vector
  ConeVertexMode vtxMode;
  TVector3 vtx_nom;               // cm
  double cosCone;                 // cos(42 deg)

  // histograms (fractions + absolute inside/outside)
  TH1D *h_fNout_pion, *h_fNout_nopion;
  TH1D *h_fQout_pion, *h_fQout_nopion;
  TH1D *h_Nin_pion, *h_Nout_pion, *h_Nin_nopion, *h_Nout_nopion;
  TH1D *h_Qin_pion, *h_Qout_pion, *h_Qin_nopion, *h_Qout_nopion;
  // same-side
  TH2D *h2_Qin_Nin_pion,  *h2_Qin_Nin_nopion;
  TH2D *h2_Qout_Nout_pion, *h2_Qout_Nout_nopion;

  // cross
  TH2D *h2_Qin_Nout_pion,  *h2_Qin_Nout_nopion;
  TH2D *h2_Qout_Nin_pion,  *h2_Qout_Nin_nopion;

  Cone42Study(const GeoCache& g, const TVector3& beam_unit, ConeVertexMode mode);

  void Fill(const EventSummary& e) override {};
  void FillEvent(WCSimRootTrigger* ev,
                 const EventSummary& e,
                 double triggerShift,
                 double triggerTime) override;

  void WritePlots() override;

private:
  TVector3 PickVertex(const EventSummary& e) const;
};

