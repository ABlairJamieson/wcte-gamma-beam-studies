#include "Studies_QvsN_2D.h"
#include "TCanvas.h"
#include "TLegend.h"

QvsN2DStudy::QvsN2DStudy() {
  h_pion   = new TH2D("h2_QvN_pion","TotQ vs NDigi (pion);NDigi;TotQ (p.e.)",70,0,1400,50,0,5000);
  h_nopion = new TH2D("h2_QvN_nopion","TotQ vs NDigi (no pion);NDigi;TotQ (p.e.)",70,0,1400,50,0,5000);
}

void QvsN2DStudy::Fill(const EventSummary& e) {
  if (e.truth.hasPions) h_pion->Fill(e.reco.nDigiHits, e.reco.totQ);
  else                 h_nopion->Fill(e.reco.nDigiHits, e.reco.totQ);
}

void QvsN2DStudy::WritePlots() {
  TCanvas c("c","c",900,700);

  c.Clear(); h_nopion->Draw("colz"); c.SaveAs("Q_vs_NHits_nopion.pdf");
  c.Clear(); h_pion->Draw("colz");   c.SaveAs("Q_vs_NHits_pion.pdf");

  TCanvas co("co","overlay",900,700);
  h_nopion->SetLineColor(kBlue);
  h_pion->SetLineColor(kRed);
  h_nopion->Draw("box");
  h_pion->Draw("box same");
  TLegend leg(0.60,0.75,0.88,0.88);
  leg.SetBorderSize(0);
  leg.AddEntry(h_nopion,"no pion","f");
  leg.AddEntry(h_pion,"pion","f");
  leg.Draw();
  co.SaveAs("Q_vs_NHits_overlay_box.pdf");
}