#pragma once
#include "Study.h"
#include "TH1D.h"

struct BasicSpectraStudy : IStudy {
  TH1D *hN_pion, *hN_nopion, *hQ_pion, *hQ_nopion, *hN_em, *hQ_em;
  TH1D *ht_pion, *ht_nopion, *ht_em;
  BasicSpectraStudy();
  void Fill(const EventSummary& e) override;
  void WritePlots() override;
};