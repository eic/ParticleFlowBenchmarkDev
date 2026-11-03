// ----------------------------------------------------------------------------
// Script for plotting results from DoubleParticleBenchmark.cxx
// (charged + neutral, Expected / Remnant clusters).
// ALL plots are drawn through DrawOverlay (TH1 / TProfile / TEfficiency).
//
// Run:     root -l -b -q Plotter.cxx
// Input:   Plots/<label>/<label>.root
// Output:  Plots/<label>/All_<label>.pdf          (all plots in one file)
//          Plots/<label>/Single/<canvas>.pdf      (one PDF per plot)
// ----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TEfficiency.h>

using namespace std;

// ============================== SETTINGS ====================================
struct Sample {
  string label;
  string title;
  double eta[2];
  string chg;
  string neu;
};

const vector<Sample> samples = {
  {"ElectronPhotonBarrel",  "e^{-}+#gamma Barrel ECal",   {-1.8,  1.5}, "e",   "#gamma"},
  {"ElectronPhotonNEndCap", "e^{-}+#gamma N Endcap ECal", {-3.6, -1.7}, "e",   "#gamma"},
  {"ElectronPhotonPEndCap", "e^{-}+#gamma P Endcap ECal", { 1.3,  3.6}, "e",   "#gamma"},
  {"PionNeutronBarrel",     "#pi^{+}+n Barrel HCal",      {-1.2,  1.2}, "#pi", "n"},
  {"PionNeutronNEndCap",    "#pi^{+}+n N Endcap HCal",    {-3.5, -1.0}, "#pi", "n"},
  {"PionNeutronPEndCap",    "#pi^{+}+n P Endcap HCal",    { 3.1,  3.8}, "#pi", "n"},
  {"PionNeutronLFHCAL",     "#pi^{+}+n LFHCAL",           { 1.0,  3.5}, "#pi", "n"},
};

const int    kColCharged = kAzure + 1;
const int    kColNeutral = kRed + 1;
const int    kColTot     = kGreen + 2;
const int    kColGray    = kGray + 2;
const double kTextSize   = 0.045;
const bool   kHeader     = true;
const int    kCanvasW    = 800;
const int    kCanvasH    = 700;
const double kDRMax      = 0.5;    // = drMax in the macro
const double kEMax       = 20.0;   // = range of *_vs_E in the macro
// ============================================================================


// ---------- helpers ----------------------------------------------------------
vector<TCanvas*> gCanvases;

// Zapisuje kazdy canvas jako osobny PDF do singleDir
// oraz wszystkie canvasy razem do jednego PDF-a (file).
void FlushCombined(const string& file, const string& singleDir) {
  if (gCanvases.empty()) return;

  // --- pojedyncze PDF-y ---
  gSystem->mkdir(singleDir.c_str(), true);
  for (auto c : gCanvases) {
    c->Print((singleDir + "/" + c->GetName() + ".pdf").c_str());
  }

  // --- jeden wspolny PDF ---
  gCanvases.front()->Print((file + "[").c_str());
  for (auto c : gCanvases) c->Print(file.c_str());
  gCanvases.back()->Print((file + "]").c_str());

  for (auto c : gCanvases) delete c;
  gCanvases.clear();
}

void SetStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillStyle(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.3, "Y");
}

template <typename T>
T* Get(TFile* f, const char* name) {
  TObject* o = f->Get(name);
  if (!o) cerr << "BRAK obiektu: " << name << " w " << f->GetName() << endl;
  return dynamic_cast<T*>(o);
}

TCanvas* NewCanvas(const Sample& s, const string& name) {
  auto c = new TCanvas((name + "_" + s.label).c_str(), "", kCanvasW, kCanvasH);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetTopMargin(kHeader ? 0.07 : 0.03);
  c->SetRightMargin(0.05);
  return c;
}

void Save(TCanvas* c) { gCanvases.push_back(c); }

void Header(const string& title) {
  if (!kHeader) return;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.045);
  t.DrawLatex(gPad->GetLeftMargin(), 0.94, ("ePIC Benchmark | " + title).c_str());
}

void RefLine(double x1, double y1, double x2, double y2) {
  auto l = new TLine(x1, y1, x2, y2);
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 2);
  l->Draw();
}


// ---------- DrawOverlay: TH1 / TProfile (b may be nullptr) -------------------
// vref >= 0: dashed vertical line (TH1) or horizontal line (TProfile).
void DrawOverlay(const Sample& s, TH1* a, const string& la, TH1* b, const string& lb,
                 const string& xt, const string& yt, const string& name,
                 double vref = -1, int colA = kColCharged, int colB = kColNeutral) {
  if (!a) return;
  const bool prof = dynamic_cast<TProfile*>(a) != nullptr;
  auto c = NewCanvas(s, name);

  vector<TH1*> hs = {a};
  vector<int>  cs = {colA};
  vector<string> ls = {la};
  if (b) { hs.push_back(b); cs.push_back(colB); ls.push_back(lb); }

  double ymax = 0;
  for (auto h : hs) ymax = std::max(ymax, h->GetMaximum());

  auto leg = new TLegend(0.2, b ? 0.78 : 0.84, 0.94, 0.90);
  leg->SetTextSize(kTextSize);

  for (size_t i = 0; i < hs.size(); ++i) {
    TH1* h = hs[i];
    h->SetTitle("");
    h->SetLineColor(cs[i]);
    h->SetLineWidth(2);
    if (prof) {
      h->SetMarkerColor(cs[i]);
      h->SetMarkerStyle(20 + i);
      h->SetMarkerSize(1.0);
      h->SetMinimum(0.0);
      h->SetMaximum(4.0);
    } else {
      h->SetFillColorAlpha(cs[i], 0.3);
    }
    if (i == 0) {
      h->GetXaxis()->SetTitle(xt.c_str());
      h->GetYaxis()->SetTitle(yt.c_str());
      if (!prof) h->SetMaximum(1.45 * ymax);
    }
    h->Draw(i == 0 ? (prof ? "E1" : "hist") : (prof ? "E1 same" : "hist same"));
    if (prof)
      leg->AddEntry(h, ls[i].c_str(), "lp");
    else
      leg->AddEntry(h, Form("%s (#mu=%.3f, #sigma=%.3f)", ls[i].c_str(), h->GetMean(), h->GetStdDev()), "f");
  }

  if (vref >= 0) {
    if (prof) RefLine(a->GetXaxis()->GetXmin(), vref, a->GetXaxis()->GetXmax(), vref);
    else      RefLine(vref, 0, vref, a->GetMaximum());
  }
  Header(s.title);
  leg->Draw();
  Save(c);
}

// ---------- DrawOverlay: TEfficiency (b may be nullptr) ----------------------
void DrawOverlay(const Sample& s, TEfficiency* a, const string& la, TEfficiency* b, const string& lb,
                 const string& xt, const string& name, double xlo, double xhi) {
  if (!a) return;
  auto c = NewCanvas(s, name);
  auto fr = c->DrawFrame(xlo, 0.0, xhi, 1.25);
  fr->GetXaxis()->SetTitle(xt.c_str());
  fr->GetYaxis()->SetTitle("Efficiency");
  auto leg = new TLegend(0.5, b ? 0.78 : 0.84, 0.94, 0.90);
  leg->SetTextSize(kTextSize);

  vector<TEfficiency*> es = {a};
  vector<int> cs = {kColCharged};
  vector<string> ls = {la};
  if (b) { es.push_back(b); cs.push_back(kColNeutral); ls.push_back(lb); }

  for (size_t i = 0; i < es.size(); ++i) {
    es[i]->SetLineColor(cs[i]); es[i]->SetMarkerColor(cs[i]);
    es[i]->SetMarkerStyle(20 + i); es[i]->SetMarkerSize(1.0); es[i]->SetLineWidth(2);
    es[i]->Draw("P same");
    leg->AddEntry(es[i], ls[i].c_str(), "lp");
  }
  RefLine(xlo, 1.0, xhi, 1.0);
  leg->Draw();
  Header(s.title);
  Save(c);
}


// ---------- all plots for one sample ----------------------------------------
void DrawSample(const Sample& s, TFile* f) {
  const string dir = "Plots/" + s.label;
  gSystem->mkdir(dir.c_str(), true);
  gCanvases.clear();
  const string pdf       = dir + "/All_" + s.label + ".pdf";
  const string singleDir = dir + "/Single";

  const string C = s.chg;
  const string N = s.neu;

  // --- load ---
  auto hEMCc  = Get<TH1D>(f, "h_EMC_Charged");
  auto hEMCn  = Get<TH1D>(f, "h_EMC_Neutral");
  auto hDRMC  = Get<TH1D>(f, "h_DeltaR_MC_charged_neutral");
  auto hDRClC = Get<TH1D>(f, "h_DeltaR_MC_Charged");
  auto hDRClN = Get<TH1D>(f, "h_DeltaR_MC_Neutral");
  auto hNExp  = Get<TH1D>(f, "h_NExpected");
  auto hNRem  = Get<TH1D>(f, "h_NRemnant");
  auto hEExp  = Get<TH1D>(f, "h_EExp");
  auto hERem  = Get<TH1D>(f, "h_ERem");
  auto rC     = Get<TH1D>(f, "h_Ratio_Charged");
  auto rN     = Get<TH1D>(f, "h_Ratio_Neutral");
  auto rT     = Get<TH1D>(f, "h_Ratio_Tot");
  auto eOp    = Get<TH1D>(f, "h_ExpOverPTrack");

  auto pCdR   = Get<TProfile>(f, "p_RatioCharged_vs_dR");
  auto pNdR   = Get<TProfile>(f, "p_RatioNeutral_vs_dR");
  auto pTdR   = Get<TProfile>(f, "p_RatioTot_vs_dR");
  auto pEPdR  = Get<TProfile>(f, "p_ExpOverP_vs_dR");
  auto pCEta  = Get<TProfile>(f, "p_RatioCharged_vs_Eta");
  auto pNEta  = Get<TProfile>(f, "p_RatioNeutral_vs_Eta");
  auto pCE    = Get<TProfile>(f, "p_RatioCharged_vs_E");
  auto pNE    = Get<TProfile>(f, "p_RatioNeutral_vs_E");

  auto efCdR  = Get<TEfficiency>(f, "eff_Charged_vs_dR");
  auto efNdR  = Get<TEfficiency>(f, "eff_Neutral_vs_dR");
  auto efCE   = Get<TEfficiency>(f, "eff_Charged_vs_E");
  auto efNE   = Get<TEfficiency>(f, "eff_Neutral_vs_E");
  auto efCEta = Get<TEfficiency>(f, "eff_Charged_vs_Eta");
  auto efNEta = Get<TEfficiency>(f, "eff_Neutral_vs_Eta");

  const string dRlab  = "#DeltaR_{MC}(" + C + "," + N + ")";
  const string etaLab = "#eta_{MC}";
  const string labC = "Charged (" + C + ")";
  const string labN = "Neutral (" + N + ")";
  const string totLab = "(E_{exp}+E_{rem}) / (E_{MC," + C + "}+E_{MC," + N + "})";

  // --- MC truth ---
  DrawOverlay(s, hEMCc, labC, hEMCn, labN, "E_{MC} [GeV]", "Events", "c_mc_energy");
  DrawOverlay(s, hDRMC, dRlab, nullptr, "", dRlab, "Events", "c_mc_deltaR", -1, kColGray);

  // --- Distance cluster <-> MC ---
  DrawOverlay(s, hDRClC, "Expected - MC " + C, hDRClN, "Remnant - MC " + N,
              "#DeltaR(cluster, MC)", "Events", "c_deltaR_cluster_MC");

  // --- Multiplicities ---
  DrawOverlay(s, hNExp, "Expected clusters", hNRem, "Remnant clusters",
              "N_{clusters}", "Events", "c_n_clusters");

  // --- Reconstructed energies ---
  DrawOverlay(s, hEExp, "E_{exp} (" + C + ")", hERem, "E_{rem} (" + N + ")",
              "E_{reco} [GeV]", "Events", "c_E_reco");

  // --- Energy response ---
  DrawOverlay(s, rC, "E_{exp}/E_{MC," + C + "}", rN, "E_{rem}/E_{MC," + N + "}",
              "E_{reco} / E_{MC}", "Events", "c_ratio_charged_neutral", 1.0);
  DrawOverlay(s, rT, totLab, nullptr, "", totLab, "Events", "c_ratio_total", 1.0, kColTot);
  DrawOverlay(s, eOp, "E_{exp} / p_{Tr}", nullptr, "", "E_{exp} / p_{Tr}", "Events",
              "c_Eexp_over_pTr", 1.0, kColCharged);

  // --- Profiles vs dR ---
  DrawOverlay(s, pCdR, "Charged: E_{exp}/E_{MC," + C + "}", pNdR, "Neutral: E_{rem}/E_{MC," + N + "}",
              dRlab, "< E_{reco} / E_{MC} >", "c_profile_ratio_vs_dR", 1.0);
  DrawOverlay(s, pTdR, totLab, nullptr, "", dRlab, "< E_{reco,tot} / E_{MC,tot} >",
              "c_profile_total_vs_dR", 1.0, kColTot);
  DrawOverlay(s, pEPdR, "E_{exp} / p_{Tr}", nullptr, "", dRlab, "< E_{exp} / p_{Tr} >",
              "c_profile_Eexp_over_pTr_vs_dR", 1.0, kColCharged);

  // --- Profiles vs eta ---
  DrawOverlay(s, pCEta, "Charged: E_{exp}/E_{MC," + C + "}", pNEta, "Neutral: E_{rem}/E_{MC," + N + "}",
              etaLab, "< E_{reco} / E_{MC} >", "c_profile_ratio_vs_eta", 1.0);

  // --- Profiles vs E ---
  DrawOverlay(s, pCE, "Charged: E_{exp}/E_{MC," + C + "}", pNE, "Neutral: E_{rem}/E_{MC," + N + "}",
              "E_{MC} [GeV]", "< E_{reco} / E_{MC} >", "c_profile_ratio_vs_E", 1.0);

  // --- Efficiencies ---
  DrawOverlay(s, efCdR, "Charged (expected)", efNdR, "Neutral (remnant)",
              dRlab, "c_eff_vs_dR", 0.0, kDRMax);
  DrawOverlay(s, efCE, "Charged (expected)", efNE, "Neutral (remnant)",
              "E_{MC} [GeV]", "c_eff_vs_E", 0.0, kEMax);
  DrawOverlay(s, efCEta, "Charged (expected)", efNEta, "Neutral (remnant)",
              etaLab, "c_eff_vs_eta", s.eta[0], s.eta[1]);

  FlushCombined(pdf, singleDir);
  cout << "  -> " << pdf << endl;
  cout << "  -> " << singleDir << "/*.pdf" << endl;
}


// ---------- main ------------------------------------------------------------
void Plotter() {
  SetStyle();
  gROOT->SetBatch(kTRUE);
  for (const auto& s : samples) {
    const string path = "Plots/" + s.label + "/" + s.label + ".root";
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { cerr << "Nie moge otworzyc " << path << "\n"; continue; }
    cout << "Rysuje: " << s.label << endl;
    DrawSample(s, f);
    f->Close();
  }
}