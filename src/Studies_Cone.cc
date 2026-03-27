#include "Studies_Cone.h"
#include "Utils.h"

#include "WCSimRootEvent.hh"

#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

#include <cmath>
#include <string>

std::string ConeStudy::Tag() const {
  // integer tag is easiest for file names + ROOT object names
  int ideg = (int)std::lround(cone_deg);
  return "Cone" + std::to_string(ideg);
}

ConeStudy::ConeStudy(const GeoCache& g,
                     const TVector3& beam_unit,
                     ConeVertexMode mode,
                     double coneDeg)
  : geo(g),
    beam_hat(beam_unit.Unit()),
    vtxMode(mode),
    vtx_nom(0.0, -42.0, -154.0),
    cone_deg(coneDeg),
    cosCone(std::cos(cone_deg * TMath::DegToRad()))
{
  const std::string tag = Tag();

  // 1D fraction overlays
  h_fNout_pion   = new TH1D((tag+"_hfNout_pion").c_str(),
                           (tag+": Nout/Ntot outside cone (pion);N_{out}/N_{tot};Events").c_str(),
                           100,0,1);
  h_fNout_nopion = new TH1D((tag+"_hfNout_nopion").c_str(),
                           (tag+": Nout/Ntot outside cone (no pion);N_{out}/N_{tot};Events").c_str(),
                           100,0,1);
  h_fQout_pion   = new TH1D((tag+"_hfQout_pion").c_str(),
                           (tag+": Qout/Qtot outside cone (pion);Q_{out}/Q_{tot};Events").c_str(),
                           100,0,1);
  h_fQout_nopion = new TH1D((tag+"_hfQout_nopion").c_str(),
                           (tag+": Qout/Qtot outside cone (no pion);Q_{out}/Q_{tot};Events").c_str(),
                           100,0,1);

  // repeat histograms for fNout_emshower and fQout_emshower studies
  h_fNout_emshower_pion   = new TH1D((tag+"_hfNout_emshower_pion").c_str(),
                                   (tag+": Nout/Ntot outside cone (emshower);N_{out}/N_{tot};Events").c_str(),
                                   100,0,1);

  h_fQout_emshower_pion   = new TH1D((tag+"_hfQout_emshower_pion").c_str(),
                                   (tag+": Qout/Qtot outside cone (emshower);Q_{out}/Q_{tot};Events").c_str(),
                                   100,0,1);

  auto mk2 = [&](const char* suffix, const char* title){
    return new TH2D((tag+"_"+suffix).c_str(), title, 70,0,1400, 50,0,5000);
  };

  h2_Qin_Nin_pion   = mk2("h2_Qin_Nin_pion",   (tag+": Q_{in} vs N_{in} (pion);N_{in};Q_{in} (p.e.)").c_str());
  h2_Qin_Nin_nopion = mk2("h2_Qin_Nin_nopion", (tag+": Q_{in} vs N_{in} (no pion);N_{in};Q_{in} (p.e.)").c_str());

  h2_Qout_Nout_pion   = mk2("h2_Qout_Nout_pion",   (tag+": Q_{out} vs N_{out} (pion);N_{out};Q_{out} (p.e.)").c_str());
  h2_Qout_Nout_nopion = mk2("h2_Qout_Nout_nopion", (tag+": Q_{out} vs N_{out} (no pion);N_{out};Q_{out} (p.e.)").c_str());

  h2_Qin_Nout_pion   = mk2("h2_Qin_Nout_pion",   (tag+": Q_{in} vs N_{out} (pion);N_{out};Q_{in} (p.e.)").c_str());
  h2_Qin_Nout_nopion = mk2("h2_Qin_Nout_nopion", (tag+": Q_{in} vs N_{out} (no pion);N_{out};Q_{in} (p.e.)").c_str());

  h2_Qout_Nin_pion   = mk2("h2_Qout_Nin_pion",   (tag+": Q_{out} vs N_{in} (pion);N_{in};Q_{out} (p.e.)").c_str());
  h2_Qout_Nin_nopion = mk2("h2_Qout_Nin_nopion", (tag+": Q_{out} vs N_{in} (no pion);N_{in};Q_{out} (p.e.)").c_str());

  h2_fQ_fN_pion = new TH2D(
    (tag+"_h2_fQ_fN_pion").c_str(),
    (tag+": Q_{out}/Q_{tot} vs N_{out}/N_{tot} (pion);N_{out}/N_{tot};Q_{out}/Q_{tot}").c_str(),
    60,0,1, 60,0,1);

  h2_fQ_fN_nopion = new TH2D(
    (tag+"_h2_fQ_fN_nopion").c_str(),
    (tag+": Q_{out}/Q_{tot} vs N_{out}/N_{tot} (no pion);N_{out}/N_{tot};Q_{out}/Q_{tot}").c_str(),
    60,0,1, 60,0,1);

}

TVector3 ConeStudy::PickVertex(WCSimRootTrigger* ev, const EventSummary& e) const {
  TVector3 vtx = vtx_nom;
  if (vtxMode == ConeVertexMode::Truth) {
    double vx, vy, vz;
    if (GetTruthVertexXYZ(ev, e.truth, vx, vy, vz)) vtx = TVector3(vx, vy, vz);
  }
  return vtx;
}

void ConeStudy::FillEvent(WCSimRootTrigger* ev,
                          const EventSummary& e,
                          double /*triggerShift*/,
                          double /*triggerTime*/)
{
  TVector3 vtx = PickVertex(ev, e);

  int Nin=0, Nout=0;
  double Qin=0.0, Qout=0.0;

  const int nDigi = ev->GetNcherenkovdigihits();
  for (int i=0; i<nDigi; i++) {
    auto* dh = (WCSimRootCherenkovDigiHit*) ev->GetCherenkovDigiHits()->At(i);
    if (!dh) continue;

    int tube = dh->GetTubeId() - 1;
    if (tube < 0 || tube >= geo.nPMT) continue;

    TVector3 u = (geo.pmt_pos[tube] - vtx).Unit();
    bool inCone = (u.Dot(beam_hat) >= cosCone);

    double q = dh->GetQ();
    // accumulate only prompt hits (e.g. within XX ns of vertex time) to avoid bias from late hits from neutrons, etc.
    // this is a somewhat loose cut, but should be sufficient for a first pass; can always tighten later if needed
    double t = dh->GetT() - geo.pmt_tof[tube]; // time minus TOF from vertex to PMT
    if (t > 100.0) continue;
    if (inCone) { Nin++;  Qin += q; }
    else        { Nout++; Qout += q; }
  }

  const int Ntot = Nin + Nout;
  const double Qtot = Qin + Qout;
  if (Ntot <= 0 || Qtot <= 0) return;

  double fN = (double)Nout / Ntot;
  double fQ = Qout / Qtot;

  if (e.truth.hasPions) {
    h_fNout_pion->Fill(fN);
    h_fQout_pion->Fill(fQ);
    h2_fQ_fN_pion->Fill(fN, fQ);
    h2_Qin_Nin_pion->Fill(Nin, Qin);
    h2_Qout_Nout_pion->Fill(Nout, Qout);
    h2_Qin_Nout_pion->Fill(Nout, Qin);
    h2_Qout_Nin_pion->Fill(Nin, Qout);
  } else {
    h_fNout_nopion->Fill(fN);
    h_fQout_nopion->Fill(fQ);
    h2_fQ_fN_nopion->Fill(fN, fQ);
    h2_Qin_Nin_nopion->Fill(Nin, Qin);
    h2_Qout_Nout_nopion->Fill(Nout, Qout);
    h2_Qin_Nout_nopion->Fill(Nout, Qin);
    h2_Qout_Nin_nopion->Fill(Nin, Qout);

    if (!e.truth.hasHadronic) {
      h_fNout_emshower_pion->Fill(fN);
      h_fQout_emshower_pion->Fill(fQ);
    }
  }
}

void ConeStudy::WritePlots()
{
  const std::string tag = Tag();

  auto overlay1D = [&](TH1* pion, TH1* nopion, const std::string& out, const char* title){
    TCanvas c("c","c",800,600);
    nopion->SetTitle(title);
    nopion->SetLineColor(kBlue); nopion->SetLineWidth(2);
    pion->SetLineColor(kRed);    pion->SetLineWidth(2);
    nopion->Draw("hist");
    pion->Draw("hist same");
    if (auto leg = c.BuildLegend(0.60,0.75,0.88,0.88)) leg->SetBorderSize(0);
    c.SaveAs(out.c_str());
  };

  auto overlay1Dem = [&](TH1* pion, TH1* nopion, TH1* emshower, const std::string& out, const char* title){
    TCanvas c("c","c",800,600);
    nopion->SetTitle(title);
    nopion->SetLineColor(kBlue); nopion->SetLineWidth(2);
    pion->SetLineColor(kRed);    pion->SetLineWidth(2);
    emshower->SetLineColor(kGreen+2); emshower->SetLineWidth(2);
    nopion->Draw("hist");
    pion->Draw("hist same");
    emshower->Draw("hist same");
    if (auto leg = c.BuildLegend(0.60,0.75,0.88,0.88)) leg->SetBorderSize(0);
    c.SaveAs(out.c_str());
  };

  auto BoxOverlay2D = [&](TH2* hNo, TH2* hPi,
                          const std::string& out, const char* title){
    TCanvas c("cbox","cbox",900,700);
    hNo->SetTitle(title);

    hNo->SetLineColor(kBlue);
    hPi->SetLineColor(kRed);

    hNo->SetLineWidth(2);
    hNo->Draw("box");
    hPi->SetLineWidth(2);
    hPi->Draw("box same");

    TLegend leg(0.60,0.75,0.88,0.88);
    leg.SetBorderSize(0);
    leg.AddEntry(hNo,"no pion","l");
    leg.AddEntry(hPi,"pion","l");
    leg.Draw();

    c.SaveAs(out.c_str());
  };

  auto mkNorm1D = [](TH1D* h, const char* name){
    TH1D* c = (TH1D*)h->Clone(name);
    double I = c->Integral(0, c->GetNbinsX()+1); // include under/overflow
    if (I > 0) c->Scale(1.0/I);
    return c;
  };

  auto Normalize2D = [](TH2* h){
    double I = h->Integral(0, h->GetNbinsX()+1, 0, h->GetNbinsY()+1);
    if (I > 0) h->Scale(1.0/I);
  };

  // 1D overlays
  overlay1D(h_fNout_pion, h_fNout_nopion,
            tag+"_fNout_overlay.pdf",
            (tag+": Outside-cone hit fraction;N_{out}/N_{tot};Events").c_str());
  overlay1D(h_fQout_pion, h_fQout_nopion,
            tag+"_fQout_overlay.pdf",
            (tag+": Outside-cone charge fraction;Q_{out}/Q_{tot};Events").c_str());

  overlay1Dem(h_fNout_pion, h_fNout_nopion, h_fNout_emshower_pion,
              tag+"_fNout_overlay_emshower.pdf",
              (tag+": Outside-cone hit fraction;N_{out}/N_{tot};Events").c_str());
  overlay1Dem(h_fQout_pion, h_fQout_nopion, h_fQout_emshower_pion,
              tag+"_fQout_overlay_emshower.pdf",
              (tag+": Outside-cone charge fraction;Q_{out}/Q_{tot};Events").c_str());

  TH1D* fNout_no_norm = mkNorm1D(h_fNout_nopion, "h_fNout_nopion_norm");
  TH1D* fNout_pi_norm = mkNorm1D(h_fNout_pion,   "h_fNout_pion_norm");

  TH1D* fQout_no_norm = mkNorm1D(h_fQout_nopion, "h_fQout_nopion_norm");
  TH1D* fQout_pi_norm = mkNorm1D(h_fQout_pion,   "h_fQout_pion_norm");

  overlay1D(fNout_pi_norm, fNout_no_norm,
          tag+"_fNout_overlay_norm.pdf",
          (tag+": Normalized outside-cone hit fraction;N_{out}/N_{tot};Normalized events").c_str());

  overlay1D(fQout_pi_norm, fQout_no_norm,
          tag+"_fQout_overlay_norm.pdf",
          (tag+": Normalized outside-cone charge fraction;Q_{out}/Q_{tot};Normalized events").c_str());

  // raw 2D overlays
  BoxOverlay2D(h2_Qin_Nin_nopion,  h2_Qin_Nin_pion,
              tag+"_Qin_vs_Nin_overlay_box.pdf",
              (tag+": Q_{in} vs N_{in};N_{in};Q_{in} (p.e.)").c_str());
  BoxOverlay2D(h2_Qout_Nout_nopion, h2_Qout_Nout_pion,
              tag+"_Qout_vs_Nout_overlay_box.pdf",
              (tag+": Q_{out} vs N_{out};N_{out};Q_{out} (p.e.)").c_str());
  BoxOverlay2D(h2_Qin_Nout_nopion, h2_Qin_Nout_pion,
              tag+"_Qin_vs_Nout_overlay_box.pdf",
              (tag+": Q_{in} vs N_{out};N_{out};Q_{in} (p.e.)").c_str());
  BoxOverlay2D(h2_Qout_Nin_nopion, h2_Qout_Nin_pion,
              tag+"_Qout_vs_Nin_overlay_box.pdf",
              (tag+": Q_{out} vs N_{in};N_{in};Q_{out} (p.e.)").c_str());

  BoxOverlay2D(h2_fQ_fN_nopion, h2_fQ_fN_pion,
             tag+"_fQ_vs_fN_overlay_box.pdf",
             (tag+": Q_{out}/Q_{tot} vs N_{out}/N_{tot};N_{out}/N_{tot};Q_{out}/Q_{tot}").c_str());

  // normalized 2D overlays
  auto mkNorm = [&](TH2D* h, const char* name){
    auto* c = (TH2D*)h->Clone((tag+"_"+name).c_str());
    Normalize2D(c);
    return c;
  };

  TH2D* QinNin_no_norm  = mkNorm(h2_Qin_Nin_nopion,  "h2_Qin_Nin_nopion_norm");
  TH2D* QinNin_pi_norm  = mkNorm(h2_Qin_Nin_pion,    "h2_Qin_Nin_pion_norm");
  TH2D* QoutNout_no_norm= mkNorm(h2_Qout_Nout_nopion,"h2_Qout_Nout_nopion_norm");
  TH2D* QoutNout_pi_norm= mkNorm(h2_Qout_Nout_pion,  "h2_Qout_Nout_pion_norm");
  TH2D* QinNout_no_norm = mkNorm(h2_Qin_Nout_nopion, "h2_Qin_Nout_nopion_norm");
  TH2D* QinNout_pi_norm = mkNorm(h2_Qin_Nout_pion,   "h2_Qin_Nout_pion_norm");
  TH2D* QoutNin_no_norm = mkNorm(h2_Qout_Nin_nopion, "h2_Qout_Nin_nopion_norm");
  TH2D* QoutNin_pi_norm = mkNorm(h2_Qout_Nin_pion,   "h2_Qout_Nin_pion_norm");

  TH2D* fQfN_no_norm = mkNorm(h2_fQ_fN_nopion, "h2_fQ_fN_nopion_norm");
  TH2D* fQfN_pi_norm = mkNorm(h2_fQ_fN_pion,   "h2_fQ_fN_pion_norm");

  BoxOverlay2D(QinNin_no_norm,  QinNin_pi_norm,
              tag+"_Qin_vs_Nin_overlay_box_norm.pdf",
              (tag+": Normalized Q_{in} vs N_{in};N_{in};Q_{in} (p.e.)").c_str());
  BoxOverlay2D(QoutNout_no_norm, QoutNout_pi_norm,
              tag+"_Qout_vs_Nout_overlay_box_norm.pdf",
              (tag+": Normalized Q_{out} vs N_{out};N_{out};Q_{out} (p.e.)").c_str());
  BoxOverlay2D(QinNout_no_norm, QinNout_pi_norm,
              tag+"_Qin_vs_Nout_overlay_box_norm.pdf",
              (tag+": Normalized Q_{in} vs N_{out};N_{out};Q_{in} (p.e.)").c_str());
  BoxOverlay2D(QoutNin_no_norm, QoutNin_pi_norm,
              tag+"_Qout_vs_Nin_overlay_box_norm.pdf",
              (tag+": Normalized Q_{out} vs N_{in};N_{in};Q_{out} (p.e.)").c_str());

  BoxOverlay2D(fQfN_no_norm, fQfN_pi_norm,
             tag+"_fQ_vs_fN_overlay_box_norm.pdf",
             (tag+": Normalized Q_{out}/Q_{tot} vs N_{out}/N_{tot};N_{out}/N_{tot};Q_{out}/Q_{tot}").c_str());


}
