#pragma once
#include "Study.h"
#include "TH2D.h"

struct QvsN2DStudy : IStudy {
  TH2D *h_pion, *h_nopion;
  QvsN2DStudy();
  void Fill(const EventSummary& e) override;
  void WritePlots() override;
};