// ----------------------------------------------------------------------------
// Simple script for plotting results from BenchmarkTrackCluster.cxx.
// Each plot is saved as a separate PDF (without TCanvas::Divide), ready for a paper.
// In addition, one combined, multi-page PDF is created for each calorimeter,
// containing all of that calorimeter's plots.
//
// Run:           root -l -b -q Plotting.cxx
// Output:        Plots/<label>/<name>_<label>.pdf    (individual plots)
//                Plots/<label>/All_<label>.pdf      (all plots for a calorimeter)
//
// All settings that need editing are in the "SETTINGS" section below.
// ----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

using namespace std;

// ============================== SETTINGS ====================================
struct Sample {
  string label;       // name (file: Plots/<label>.root, directory: Plots/<label>/)
  string title;       // title displayed on plots
  int    color;       // histogram color
  double eta[2];      // eta range for plots
  double pt[2];       // pT (and p) range for plots
};

const vector<Sample> samples = {
  // label            title                       color       eta            pT
  {"ElectronBarrel",  "e^{-} Barrel ECal",        kAzure+1,  {-1.8,  1.5},  {0.0, 20.0}},
  {"ElectronNEndCap", "e^{-} N Endcap ECal",      kAzure+1,  {-3.6, -1.9},  {0.0, 7.0}},
  {"ElectronPEndCap", "e^{-} P Endcap ECal",      kAzure+1,  { 1.3,  3.6},  {0.0, 10.0}},
  {"PionBarrel",      "#pi^{-} Barrel HCal",      kRed+1,    {-1.2,  1.2},  {0.0, 20.0}},
  {"PionNEndCap",     "#pi^{-} N Endcap HCal",    kRed+1,    {-3.5, -1.0},  {0.0, 13.0}},
  {"PionPEndCap",     "#pi^{-} P Endcap Insert",  kRed+1,    { 3.1,  3.8},  {0.0,  2.0}},
  {"PionLFHCAL",      "#pi^{-} LFHCAL",           kRed+1,    { 1.0,  3.5},  {0.0, 13.0}},
};

const double kTextSize  = 0.045;  // TLatex text size
const bool   kHeader    = true;   // "ePIC Benchmark | ..." header (set false for a paper)
const int    kCanvasW   = 800;    // single-plot canvas size
const int    kCanvasH   = 700;
const bool   kCombined  = true;   // create an additional All_<label>.pdf with all calorimeter plots
// ============================================================================


// ---------- helper functions -------------------------------------------------
string gSamplePdf;               // multi-page PDF for the current calorimeter
vector<TCanvas*> gCanvases;      // canvases for the current calorimeter (for the combined PDF)

// Print all collected canvases to one multi-page PDF,
// only AFTER saving the individual PDFs (to avoid interleaving writes).
void FlushCombined(const string& file) {
  if (gCanvases.empty()) return;
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
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(64);
}

string Fmt(const char* fmt, double a, double b = 0) { return string(Form(fmt, a, b)); }

TH1* Get(TFile* f, const char* name) {
  TObject* o = f->Get(name);
  if (!o) cerr << "BRAK histogramu: " << name << " w " << f->GetName() << endl;
  return (TH1*)o;
}

// Create a new canvas with margins.
TCanvas* NewCanvas(const Sample& s, const string& name, bool is2D = false) {
  auto c = new TCanvas((name + "_" + s.label).c_str(), "", kCanvasW, kCanvasH);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetTopMargin(kHeader ? 0.07 : 0.03);
  c->SetRightMargin(is2D ? 0.15 : 0.05);
  return c;
}

// Save to PDF (a separate file and a page in All_<label>.pdf), then clean up.
void Save(TCanvas* c, const Sample& s, const string& name) {
  c->SaveAs(Form("Plots/%s/%s_%s.pdf", s.label.c_str(), name.c_str(), s.label.c_str()));
  if (kCombined) gCanvases.push_back(c);   // zostaje do FlushCombined()
  else delete c;
}

void Header(const string& title) {
  if (!kHeader) return;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.045);
  t.DrawLatex(gPad->GetLeftMargin(), 0.94, ("ePIC Benchmark | " + title).c_str());
}

// Draw text lines (TLatex) in the upper-right corner.
void Text(const vector<string>& lines, double x = 0.66, double y = 0.88, int color= kBlack, bool outline = false) {
   TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(kTextSize);
  const double d = 0.002;   // grubość obwódki (współrzędne NDC)
  for (auto& l : lines) {
    if (outline) {
      t.SetTextColor(kBlack);
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
          if (dx || dy) t.DrawLatex(x + dx * d, y + dy * d, l.c_str());
    }
    t.SetTextColor(color);
    t.DrawLatex(x, y, l.c_str());
    y -= 1.3 * kTextSize;
  }
}

void Line(TH1* h, int col) {
  h->SetTitle("");
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetFillColorAlpha(col, 0.3);
}

// Fraction of particles matched within the specified range.
double HitFraction(TH1* all, TH1* hit, double lo, double hi) {
  int b1 = all->FindBin(lo), b2 = all->FindBin(hi - 1e-6);
  double n = all->Integral(b1, b2);
  return n > 0 ? hit->Integral(b1, b2) / n : 0;
}


// ---------- plot types -------------------------------------------------------

// MC truth vs. MC matched (two histograms on one plot).
void DrawTruthHit(const Sample& s, TH1* all, TH1* hit, const char* xt,
                  double lo, double hi, const string& name, double frac) {
  auto c = NewCanvas(s, name);
  Line(all, kGray + 2);
  Line(hit, s.color);
  all->GetXaxis()->SetTitle(xt);
  all->GetYaxis()->SetTitle("Particles");
  all->GetXaxis()->SetRangeUser(lo, hi);
  all->SetMinimum(0);
  all->SetMaximum(1.45 * all->GetMaximum());
  all->Draw("hist");
  hit->Draw("hist same");
  Header(s.title);
  auto leg = new TLegend(0.62, 0.78, 0.94, 0.90);
  leg->SetTextSize(kTextSize);
  leg->AddEntry(all, "MC truth", "f");
  leg->AddEntry(hit, "MC matched", "f");
  leg->Draw();
  Text({Fmt("Matched: %.1f %%", 100 * frac),
        Fmt("%.1f < #eta < %.1f", s.eta[0], s.eta[1])}, 0.65, 0.74);
  Save(c, s, name);
}

// Standard 1D histogram; stats = true adds Mean/Std/Entries.
void DrawHist(const Sample& s, TH1* h, const char* xt, const char* yt,
              const string& name, bool stats = true) {
  auto c = NewCanvas(s, name);
  Line(h, s.color);
  h->GetXaxis()->SetTitle(xt);
  h->GetYaxis()->SetTitle(yt);
  h->SetMaximum(1.45 * h->GetMaximum());
  h->Draw("hist");
  Header(s.title);
  if (stats)
    Text({Fmt("Mean = %.3f", h->GetMean()),
          Fmt("Std = %.3f",  h->GetStdDev()),
          Fmt("Entries = %.0f", h->GetEntries())});
  Save(c, s, name);
}

// 2D map, optionally with a profile (prof) and f_sub (f_sub > 0).
void Draw2D(const Sample& s, TH2* h, TH1* prof, const char* xt, const char* yt,
            const string& name, double xlo, double xhi, double f_sub = -1) {
  auto c = NewCanvas(s, name, true);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xt);
  h->GetYaxis()->SetTitle(yt);
  if (xhi > xlo) h->GetXaxis()->SetRangeUser(xlo, xhi);
  h->Draw("COLZ");
  if (prof) {
    prof->SetLineColor(kBlack);   prof->SetMarkerColor(kBlack);
    prof->SetMarkerStyle(20);     prof->SetMarkerSize(0.9); prof->SetLineWidth(2);
    if (xhi > xlo) prof->GetXaxis()->SetRangeUser(xlo, xhi);
    prof->Draw("same");
  }
  Header(s.title);
  if (prof && f_sub > 0) Text({"Black: mean profile", Fmt("f_{sub} = %.3f", f_sub)}, 0.18, 0.88, kWhite, true);
  Save(c, s, name);
}


// ---------- all plots for one sample ----------------------------------------
void DrawSample(const Sample& s, TFile* f) {
  gSystem->mkdir(("Plots/" + s.label).c_str(), true);

  TH1* mcP   = Get(f, "MCMomentum");   TH1* mcPh  = Get(f, "MCMomentumHited");
  TH1* mcEta = Get(f, "MCEta");        TH1* mcEh  = Get(f, "MCEtaHited");
  TH1* mcPt  = Get(f, "MCPt");         TH1* mcPth = Get(f, "MCPtHited");
  TH1* mcPhi = Get(f, "MCPhi");        TH1* mcPhh = Get(f, "MCPhiHited");
  TH2* dEP   = (TH2*)Get(f, "DeltaEtaPhiHist");
  TH1* rEP   = Get(f, "REtaPhiHist");
  TH1* eCoT  = Get(f, "EClusteroverETrack");
  TH1* eMCoC = Get(f, "EMCoverECluster");
  TH1* eMCoT = Get(f, "EMCoverETrack");
  TH2* r2Pt  = (TH2*)Get(f, "EnergyRatiovsTrackPt");
  TH2* r2Eta = (TH2*)Get(f, "EnergyRatiovsTrackEta");
  TH1* prPt  = Get(f, "EnergyRatioProfilevsTrackPt");
  TH1* prEta = Get(f, "EnergyRatioProfilevsTrackEta");
  TH1* nTrk  = Get(f, "NTracksMatchedToCluster");
  TH1* nCl   = Get(f, "NClustersMatchedToTrack");
  TH1* sub   = Get(f, "EClustMinusFsubEtrk");
  TH1* eEcHcT = (TH1*)f->Get("EEcalHcalOverETrack");

  if (!mcP || !mcPh || !mcEta || !mcEh || !mcPt || !mcPth || !mcPhi || !mcPhh || !dEP || !rEP ||
      !eCoT || !eMCoC || !eMCoT || !r2Pt || !r2Eta || !prPt || !prEta || !nTrk || !nCl || !sub) {
    cerr << "Pomijam " << s.label << " (brakuje histogramow)\n";
    return;
  }

  gSamplePdf = "Plots/" + s.label + "/All_" + s.label + ".pdf";
  gCanvases.clear();

  const double f_sub = eCoT->GetMean();   // same value as f_sub in the macro
  const double frac  = HitFraction(mcEta, mcEh, s.eta[0], s.eta[1]);

  // MC truth vs matched
  DrawTruthHit(s, mcP,   mcPh,  "p [GeV]",     0,  20,  "c_mc_p",   frac);
  DrawTruthHit(s, mcEta, mcEh,  "#eta",        s.eta[0], s.eta[1], "c_mc_eta", frac);
  DrawTruthHit(s, mcPt,  mcPth, "p_{T} [GeV]", s.pt[0],  s.pt[1],  "c_mc_pt",  frac);
  DrawTruthHit(s, mcPhi, mcPhh, "#phi [rad]",  -3.15,    3.15,     "c_mc_phi", frac);

  // Track-cluster matching
  Draw2D(s, dEP, nullptr, "#Delta#eta", "#Delta#phi [rad]", "c_delta_eta_phi", 0, 0);
  DrawHist(s, rEP, "R_{#eta#phi} = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}}", "Clusters", "c_R_eta_phi", false);

  // Energy ratios
  DrawHist(s, eCoT,  "E_{Cl} / p_{Tr}", "Events", "c_ECl_over_pTr");
  DrawHist(s, eMCoC, "E_{MC} / E_{Cl}", "Events", "c_EMC_over_ECl");
  DrawHist(s, eMCoT, "E_{MC} / p_{Tr}", "Events", "c_EMC_over_pTr");

  // E/p vs. pT and eta (with profiles)
  Draw2D(s, r2Pt,  prPt,  "Track p_{T} [GeV]", "E_{Cl} / E_{Tr}", "c_ratio_vs_pt",  s.pt[0],  s.pt[1],  f_sub);
  Draw2D(s, r2Eta, prEta, "Track #eta",        "E_{Cl} / E_{Tr}", "c_ratio_vs_eta", s.eta[0], s.eta[1], f_sub);

  // Multiplicities
  DrawHist(s, nTrk, "N_{tracks} per cluster", "Events", "c_ntracks_per_cluster", false);
  DrawHist(s, nCl,  "N_{clusters} per track", "Events", "c_nclusters_per_track", false);

  // E_clust - f_sub * E_trk
  {
    auto c = NewCanvas(s, "c_EClustMinusFsubEtrk");
    Line(sub, s.color);
    sub->GetXaxis()->SetTitle("E_{clust} - f_{sub} E_{trk} [GeV]");
    sub->GetYaxis()->SetTitle("Events");
    sub->SetMaximum(1.45 * sub->GetMaximum());
    sub->Draw("hist");
    Header(s.title);
    Text({Fmt("f_{sub} = %.3f", f_sub),
          Fmt("Mean = %.3f GeV", sub->GetMean()),
          Fmt("Std = %.3f GeV",  sub->GetStdDev()),
          Fmt("Entries = %.0f",  sub->GetEntries())},0.60, 0.88);
    Save(c, s, "c_EClustMinusFsubEtrk");
  }
  // (E_ECal + E_HCal) / p_Tr - pions only
  if (eEcHcT) DrawHist(s, eEcHcT, "(E_{ECal} + E_{HCal}) / p_{Tr}", "Events", "c_EEcalHcal_over_pTr");

  // Combined PDF for this calorimeter (all plots together).
  if (kCombined) FlushCombined(gSamplePdf);

}


// ---------- main ------------------------------------------------------------
void Plotting() {
  SetStyle();
  gROOT->SetBatch(kTRUE);
  for (const auto& s : samples) {
    TFile* f = TFile::Open(("Plots/" + s.label + ".root").c_str(), "READ");
    if (!f || f->IsZombie()) { cerr << "Nie moge otworzyc Plots/" << s.label << ".root\n"; continue; }
    cout << "Rysuje: " << s.label << endl;
    DrawSample(s, f);
    f->Close();
  }
}