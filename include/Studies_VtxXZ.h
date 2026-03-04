#pragma once
#include "Study.h"
#include "TH2D.h"

struct VtxXZStudy : IStudy {
  TH2D *h_all, *h_pion, *h_nopion;
  VtxXZStudy();
  void Fill(const EventSummary& e) override;
  void WritePlots() override;
};