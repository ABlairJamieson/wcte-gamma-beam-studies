#pragma once
#include "Types.h"

class WCSimRootTrigger; 

struct IStudy {
  virtual ~IStudy() = default;

  // Always available: cheap per-event summary
  virtual void Fill(const EventSummary& e) = 0;

  // Optional hook: raw event (default no-op)
  virtual void FillEvent(WCSimRootTrigger* /*ev*/,
                         const EventSummary& /*e*/,
                         double /*triggerShift*/,
                         double /*triggerTime*/) {}
  virtual void WritePlots() = 0;
};