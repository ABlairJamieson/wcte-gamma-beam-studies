#pragma once
#include <vector>
#include "TVector3.h"


struct GeoCache {
  int nPMT = 0;
  double vg_cm_per_ns = 0.0;
  double max_r = 0.0;
  double max_z = 0.0;
  double barrelCut = 0.0;

  std::vector<TVector3> pmt_pos;
  std::vector<double>   pmt_tof;
  std::vector<int>      pmt_cylLoc;
  std::vector<int>      mpmtNo;
  std::vector<int>      mpmt_pmtNo;
  std::vector<std::vector<double>> eventXY; // unrolled display coords (x,y)
};

struct EventTruth {
  bool hasHadronic = false;
  bool hasPions = false;
  bool hasPi0 = false;
  bool hasPiCh = false;
};

struct EventReco {
  int    nDigiHits = 0; // total number of digi hits in the event
  double totQ = 0.0;    // total charge in p.e. summed over all digi hits
  double t = 0.0;       // time of first digi hit (relative to trigger time, in ns)
  double tavg = 0.0;     // average time of digi hits (relative to trigger time, in ns)
};

struct EventSummary {
  long iev = -1;

  EventTruth truth;
  EventReco  reco;

  double qPerHit = 0.0;

  // truth vertex proxy (xz)
  bool hasVtxXZ = false;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  // beam axis (unit) is global; vertex mode selects cone vertex
  bool hasTruthVtx3 = false;

  int Nin=0, Nout=0;
  double Qin=0.0, Qout=0.0;
  double fNout=0.0, fQout=0.0;

};