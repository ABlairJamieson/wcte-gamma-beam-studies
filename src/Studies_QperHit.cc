#include "Studies_QperHit.h"
#include "TCanvas.h"
#include "TLegend.h"

QperHitStudy::QperHitStudy() {
  h_pion   = new TH1D("hQperHit_pion","Q/NHits (pion);Q/NHits (p.e.);Events",200,0,10);
  h_nopion = new TH1D("hQperHit_nopion","Q/NHits (no pion);Q/NHits (p.e.);Events",200,0,10);
}

void QperHitStudy::Fill(const EventSummary& e) {
  if (e.reco.nDigiHits <= 0) return;
  if (e.truth.hasPions) h_pion->Fill(e.qPerHit);
  else                 h_nopion->Fill(e.qPerHit);
}

void QperHitStudy::WritePlots() {
  TCanvas c("c","c",800,600);
  h_nopion->SetLineColor(kBlue); h_nopion->SetLineWidth(2);
  h_pion->SetLineColor(kRed);    h_pion->SetLineWidth(2);
  h_nopion->Draw("hist");
  h_pion->Draw("hist same");
  c.BuildLegend(0.60,0.75,0.88,0.88)->SetBorderSize(0);
  c.SaveAs("QperHit_overlay_raw.pdf");
}