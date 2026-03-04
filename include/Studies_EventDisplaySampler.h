#pragma once
#include "Study.h"
#include "Types.h"

struct EventDisplaySamplerStudy : IStudy {
  const GeoCache& geo;
  int nSaveEach;
  int savedGamma=0, savedPion=0;

  // we need per-event PMT charges, so pass them via EventSummary? easiest: regenerate in Run.cc
  // We'll implement this as a no-op here and handle sampling in Run.cc for now.

  EventDisplaySamplerStudy(const GeoCache& g, int n=5) : geo(g), nSaveEach(n) {}
  void Fill(const EventSummary&) override {}
  void WritePlots() override {}
};