#include "Studies_VtxXZ.h"
#include "TCanvas.h"

VtxXZStudy::VtxXZStudy() {
  h_all    = new TH2D("hVtxXZ_all","Truth vertex proxy XZ (all);x (cm);z (cm)",200,-500,500,200,-500,500);
  h_pion   = new TH2D("hVtxXZ_pion","Truth vertex proxy XZ (pion);x (cm);z (cm)",200,-500,500,200,-500,500);
  h_nopion = new TH2D("hVtxXZ_nopion","Truth vertex proxy XZ (no pion);x (cm);z (cm)",200,-500,500,200,-500,500);
}

void VtxXZStudy::Fill(const EventSummary& e) {
  if (!e.hasVtxXZ) return;
  h_all->Fill(e.vx, e.vz);
  if (e.truth.hasPions) h_pion->Fill(e.vx, e.vz);
  else                 h_nopion->Fill(e.vx, e.vz);
}

void VtxXZStudy::WritePlots() {
  TCanvas c("c","c",900,800);
  c.Clear(); h_all->Draw("colz");    c.SaveAs("VtxXZ_all.pdf");
  c.Clear(); h_pion->Draw("colz");   c.SaveAs("VtxXZ_pion.pdf");
  c.Clear(); h_nopion->Draw("colz"); c.SaveAs("VtxXZ_nopion.pdf");
}