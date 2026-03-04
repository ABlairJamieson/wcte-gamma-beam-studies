#pragma once
#include "Study.h"
#include "TH1D.h"

struct QperHitStudy : IStudy {
  TH1D *h_pion, *h_nopion;
  QperHitStudy();
  void Fill(const EventSummary& e) override;
  void WritePlots() override;
};