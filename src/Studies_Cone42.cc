#include "Studies_Cone42.h"
#include "Utils.h"

#include "WCSimRootEvent.hh"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"


#include <cmath>

Cone42Study::Cone42Study(const GeoCache& g, const TVector3& beam_unit, ConeVertexMode mode)
  : geo(g), beam_hat(beam_unit.Unit()), vtxMode(mode),
    vtx_nom(0.0, -42.0, -154.0),
    cosCone(std::cos(42.0 * TMath::DegToRad()))
{
  h_fNout_pion   = new TH1D("hfNout_pion",   "Nout/Ntot outside 42deg (pion);N_{out}/N_{tot};Events",100,0,1);
  h_fNout_nopion = new TH1D("hfNout_nopion", "Nout/Ntot outside 42deg (no pion);N_{out}/N_{tot};Events",100,0,1);
  h_fQout_pion   = new TH1D("hfQout_pion",   "Qout/Qtot outside 42deg (pion);Q_{out}/Q_{tot};Events",100,0,1);
  h_fQout_nopion = new TH1D("hfQout_nopion", "Qout/Qtot outside 42deg (no pion);Q_{out}/Q_{tot};Events",100,0,1);

  auto mk2 = [](const char* n, const char* t){
    return new TH2D(n,t, 70,0,1400, 50,0,5000);
  };

  h2_Qin_Nin_pion   = mk2("h2_Qin_Nin_pion",   "Q_{in} vs N_{in} (pion);N_{in};Q_{in} (p.e.)");
  h2_Qin_Nin_nopion = mk2("h2_Qin_Nin_nopion", "Q_{in} vs N_{in} (no pion);N_{in};Q_{in} (p.e.)");

  h2_Qout_Nout_pion   = mk2("h2_Qout_Nout_pion",   "Q_{out} vs N_{out} (pion);N_{out};Q_{out} (p.e.)");
  h2_Qout_Nout_nopion = mk2("h2_Qout_Nout_nopion", "Q_{out} vs N_{out} (no pion);N_{out};Q_{out} (p.e.)");

  // cross
  h2_Qin_Nout_pion   = mk2("h2_Qin_Nout_pion",   "Q_{in} vs N_{out} (pion);N_{out};Q_{in} (p.e.)");
  h2_Qin_Nout_nopion = mk2("h2_Qin_Nout_nopion", "Q_{in} vs N_{out} (no pion);N_{out};Q_{in} (p.e.)");

  h2_Qout_Nin_pion   = mk2("h2_Qout_Nin_pion",   "Q_{out} vs N_{in} (pion);N_{in};Q_{out} (p.e.)");
  h2_Qout_Nin_nopion = mk2("h2_Qout_Nin_nopion", "Q_{out} vs N_{in} (no pion);N_{in};Q_{out} (p.e.)");
}

void Cone42Study::FillEvent(WCSimRootTrigger* ev,
                            const EventSummary& e,
                            double /*triggerShift*/,
                            double /*triggerTime*/)
{
  // choose vertex
  TVector3 vtx = vtx_nom;
  if (vtxMode == ConeVertexMode::Truth) {
    double vx, vy, vz;
    if (GetTruthVertexXYZ(ev, e.truth, vx, vy, vz)) vtx = TVector3(vx, vy, vz);
  }

  int Nin=0, Nout=0;
  double Qin=0.0, Qout=0.0;

  const int nDigi = ev->GetNcherenkovdigihits();
  for (int i=0; i<nDigi; i++) {
    auto* dh = (WCSimRootCherenkovDigiHit*) ev->GetCherenkovDigiHits()->At(i);
    if (!dh) continue;

    int tube = dh->GetTubeId() - 1;
    if (tube < 0 || tube >= geo.nPMT) continue;

    // direction from vertex to PMT
    TVector3 u = (geo.pmt_pos[tube] - vtx).Unit();
    bool inCone = (u.Dot(beam_hat) >= cosCone);

    double q = dh->GetQ();
    if (inCone) { Nin++;  Qin += q; }
    else        { Nout++; Qout += q; }
  }

  const int Ntot = Nin + Nout;
  const double Qtot = Qin + Qout;
  if (Ntot <= 0 || Qtot <= 0) return;

  if (e.truth.hasPions) {
    h_fNout_pion->Fill((double)Nout / Ntot);
    h_fQout_pion->Fill(Qout / Qtot);
      h2_Qin_Nin_pion->Fill(Nin, Qin);
      h2_Qout_Nout_pion->Fill(Nout, Qout);
      h2_Qin_Nout_pion->Fill(Nout, Qin);
      h2_Qout_Nin_pion->Fill(Nin, Qout);
  } else {
    h_fNout_nopion->Fill((double)Nout / Ntot);
    h_fQout_nopion->Fill(Qout / Qtot);
      h2_Qin_Nin_nopion->Fill(Nin, Qin);
      h2_Qout_Nout_nopion->Fill(Nout, Qout);
      h2_Qin_Nout_nopion->Fill(Nout, Qin);
      h2_Qout_Nin_nopion->Fill(Nin, Qout);
  }
}

void Cone42Study::WritePlots()
{
  auto overlay = [](TH1* pion, TH1* nopion, const char* out, const char* title){
    TCanvas c("c","c",800,600);
    nopion->SetTitle(title);
    nopion->SetLineColor(kBlue); nopion->SetLineWidth(2);
    pion->SetLineColor(kRed);    pion->SetLineWidth(2);
    nopion->Draw("hist");
    pion->Draw("hist same");
    c.BuildLegend(0.60,0.75,0.88,0.88)->SetBorderSize(0);
    c.SaveAs(out);
  };

  auto BoxOverlay2D = [](TH2* hNo, TH2* hPi,
                       const char* out, const char* title){
    TCanvas c("cbox","cbox",900,700);
    hNo->SetTitle(title);

    hNo->SetLineColor(kBlue);
    hPi->SetLineColor(kRed);
    //hNo->SetFillColor(kBlue);
    //hPi->SetFillColor(kRed);

    hNo->SetLineWidth(2);                    
    hNo->Draw("box");
    hPi->SetLineWidth(1);
    hPi->Draw("box same");

    TLegend leg(0.60,0.75,0.88,0.88);
    leg.SetBorderSize(0);
    leg.AddEntry(hNo,"no pion","l");
    leg.AddEntry(hPi,"pion","l");
    leg.Draw();

    c.SaveAs(out);
  };

  auto Normalize2D = [](TH2* h){
    double I = h->Integral(0, h->GetNbinsX()+1, 0, h->GetNbinsY()+1); // incl under/overflow
    if (I > 0) h->Scale(1.0/I);
  };

  overlay(h_fNout_pion, h_fNout_nopion, "Cone42_fNout_overlay.pdf",
          "Outside-cone hit fraction;N_{out}/N_{tot};Events");
  overlay(h_fQout_pion, h_fQout_nopion, "Cone42_fQout_overlay.pdf",
          "Outside-cone charge fraction;Q_{out}/Q_{tot};Events");

  BoxOverlay2D(h2_Qin_Nin_nopion,  h2_Qin_Nin_pion,
             "Cone42_Qin_vs_Nin_overlay_box.pdf",
             "Q_{in} vs N_{in};N_{in};Q_{in} (p.e.)");

  BoxOverlay2D(h2_Qout_Nout_nopion, h2_Qout_Nout_pion,
             "Cone42_Qout_vs_Nout_overlay_box.pdf",
             "Q_{out} vs N_{out};N_{out};Q_{out} (p.e.)");

  BoxOverlay2D(h2_Qin_Nout_nopion, h2_Qin_Nout_pion,
             "Cone42_Qin_vs_Nout_overlay_box.pdf",
             "Q_{in} vs N_{out};N_{out};Q_{in} (p.e.)");

  BoxOverlay2D(h2_Qout_Nin_nopion, h2_Qout_Nin_pion,
             "Cone42_Qout_vs_Nin_overlay_box.pdf",
             "Q_{out} vs N_{in};N_{in};Q_{out} (p.e.)");

  // Qin vs Nin
  TH2D* h2_Qin_Nin_nopion_norm = (TH2D*)h2_Qin_Nin_nopion->Clone("h2_Qin_Nin_nopion_norm");
  TH2D* h2_Qin_Nin_pion_norm   = (TH2D*)h2_Qin_Nin_pion->Clone("h2_Qin_Nin_pion_norm");
  Normalize2D(h2_Qin_Nin_nopion_norm);
  Normalize2D(h2_Qin_Nin_pion_norm);

  BoxOverlay2D(h2_Qin_Nin_nopion_norm, h2_Qin_Nin_pion_norm,
             "Cone42_Qin_vs_Nin_overlay_box_norm.pdf",
             "Normalized Q_{in} vs N_{in};N_{in};Q_{in} (p.e.)");


  // Qout vs Nout
  TH2D* h2_Qout_Nout_nopion_norm = (TH2D*)h2_Qout_Nout_nopion->Clone("h2_Qout_Nout_nopion_norm");
  TH2D* h2_Qout_Nout_pion_norm   = (TH2D*)h2_Qout_Nout_pion->Clone("h2_Qout_Nout_pion_norm");
  Normalize2D(h2_Qout_Nout_nopion_norm);
  Normalize2D(h2_Qout_Nout_pion_norm);

  BoxOverlay2D(h2_Qout_Nout_nopion_norm, h2_Qout_Nout_pion_norm,
             "Cone42_Qout_vs_Nout_overlay_box_norm.pdf",
             "Normalized Q_{out} vs N_{out};N_{out};Q_{out} (p.e.)");


  // Qin vs Nout
  TH2D* h2_Qin_Nout_nopion_norm = (TH2D*)h2_Qin_Nout_nopion->Clone("h2_Qin_Nout_nopion_norm");
  TH2D* h2_Qin_Nout_pion_norm   = (TH2D*)h2_Qin_Nout_pion->Clone("h2_Qin_Nout_pion_norm");
  Normalize2D(h2_Qin_Nout_nopion_norm);
  Normalize2D(h2_Qin_Nout_pion_norm);

  BoxOverlay2D(h2_Qin_Nout_nopion_norm, h2_Qin_Nout_pion_norm,
             "Cone42_Qin_vs_Nout_overlay_box_norm.pdf",
             "Normalized Q_{in} vs N_{out};N_{out};Q_{in} (p.e.)");


  // Qout vs Nin
  TH2D* h2_Qout_Nin_nopion_norm = (TH2D*)h2_Qout_Nin_nopion->Clone("h2_Qout_Nin_nopion_norm");
  TH2D* h2_Qout_Nin_pion_norm   = (TH2D*)h2_Qout_Nin_pion->Clone("h2_Qout_Nin_pion_norm");
  Normalize2D(h2_Qout_Nin_nopion_norm);
  Normalize2D(h2_Qout_Nin_pion_norm);

  BoxOverlay2D(h2_Qout_Nin_nopion_norm, h2_Qout_Nin_pion_norm,
             "Cone42_Qout_vs_Nin_overlay_box_norm.pdf",
             "Normalized Q_{out} vs N_{in};N_{in};Q_{out} (p.e.)");
}