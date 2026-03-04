#pragma once
#include <string>
#include <vector>
#include "Types.h"

class TFile;
class TTree;
class TH1;
class TH1D;
class WCSimRootTrigger;

bool GetPrimaryVtxAndDir(WCSimRootTrigger* ev, TVector3& vtx, TVector3& dir);
bool BuildGeoCache(TFile* f, const TVector3& vtx, GeoCache& geo);

EventTruth ClassifyTruth(WCSimRootTrigger* ev);

EventReco ComputeReco(WCSimRootTrigger* ev, const GeoCache& geo,
                      std::vector<double>* pmt_charge_accum,
                      TH1D* hDigiTimeMinusTOF,
                      double triggerShift, double triggerTime,
                      double timeCutNs = 1e9);

bool GetTruthVertexXZ(WCSimRootTrigger* ev, const EventTruth& truth,
                      double& vx, double& vz);

bool GetTruthVertexXYZ(WCSimRootTrigger* ev, const EventTruth& truth,
                       double& vx, double& vy, double& vz);

void NormalizeToUnitArea(TH1* h);

void SaveEventDisplayPDF(const GeoCache& geo,
                         const std::vector<double>& pmt_charge,
                         const std::string& outname);
