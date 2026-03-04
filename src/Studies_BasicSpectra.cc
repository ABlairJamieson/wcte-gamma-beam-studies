#include "Studies_BasicSpectra.h"
#include "TCanvas.h"
#include "TLegend.h"

BasicSpectraStudy::BasicSpectraStudy() {
  hN_pion   = new TH1D("hN_pion","NDigi (pion);NDigi;Events",2000,0,2000);
  hN_nopion = new TH1D("hN_nopion","NDigi (no pion);NDigi;Events",2000,0,2000);
  hQ_pion   = new TH1D("hQ_pion","TotQ (pion);TotQ (p.e.);Events",10000,0,10000);
  hQ_nopion = new TH1D("hQ_nopion","TotQ (no pion);TotQ (p.e.);Events",10000,0,10000);
}

void BasicSpectraStudy::Fill(const EventSummary& e) {
  if (e.truth.hasPions) { hN_pion->Fill(e.reco.nDigiHits); hQ_pion->Fill(e.reco.totQ); }
  else                 { hN_nopion->Fill(e.reco.nDigiHits); hQ_nopion->Fill(e.reco.totQ); }
}

static void Overlay(TH1* pion, TH1* nopion, const char* out, const char* title, bool logy=false) {
  TCanvas c("c","c",800,600);
  nopion->SetTitle(title);
  nopion->SetLineColor(kBlue); nopion->SetLineWidth(2);
  pion->SetLineColor(kRed);    pion->SetLineWidth(2);
  nopion->Draw("hist");
  pion->Draw("hist same");
  auto leg = c.BuildLegend(0.60,0.75,0.88,0.88); leg->SetBorderSize(0);
  if (logy) c.SetLogy();
  c.SaveAs(out);
}

void BasicSpectraStudy::WritePlots() {
  Overlay(hN_pion, hN_nopion, "NDigi_overlay_raw.pdf", "NDigi: pion vs no-pion;NDigi;Events", true);
  Overlay(hQ_pion, hQ_nopion, "TotQ_overlay_raw.pdf",  "TotQ: pion vs no-pion;TotQ (p.e.);Events", true);
}