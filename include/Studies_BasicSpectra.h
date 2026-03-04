#pragma once
#include "Study.h"
#include "TH1D.h"

struct BasicSpectraStudy : IStudy {
  TH1D *hN_pion, *hN_nopion, *hQ_pion, *hQ_nopion;
  BasicSpectraStudy();
  void Fill(const EventSummary& e) override;
  void WritePlots() override;
};