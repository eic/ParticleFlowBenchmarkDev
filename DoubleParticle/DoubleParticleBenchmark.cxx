// ----------------------------------------------------------------------------
// Macro for analyzing charged + neutral events using the Expected / Remnant clusters
//
//   *ExpectedClusters -> charged particle (cluster built from the track projection)
//   *RemnantClusters  -> neutral particle (what remains after the charged particle is removed)
//
// No own cluster-track matching is done here. The track is used only as a
// reference for the energy (p_track taken from the track projection that
// the Expected cluster sits on).
//
// Usage:
//   root -b -q DoubleParticleBenchmark.cxx++
// ----------------------------------------------------------------------------

#define BENCHMARKCHARGEDNEUTRAL_CXX

#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TMath.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/TrackPoint.h>
#include <edm4eic/TrackSegment.h>
#include <edm4eic/TrackSegmentCollection.h>

#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticleCollection.h>

using std::string;
using std::cout;
using std::endl;
using std::numeric_limits;

struct Options {
  string in_file;
  string out_file;
  string expectedClusters;          // charged particle   (e.g. EcalEndcapNExpectedClusters)
  string remnantClusters;           // neutral particle     (e.g. EcalEndcapNRemnantClusters)
  string projections;               // track projections (for p_track)
  uint64_t system;                  // calorimeter system ID (for projections)
  bool debug;
  std::array<double, 2> etaRange;
  int pdgExpected;                  // charged particle abs(PDG)
  int pdgRemnant;                   // neutral particle abs(PDG)
};

const Options ElectronPhotonNEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonNEndCap/ElectronPhotonNEndCap.root",
  "EcalEndcapNExpectedClusters",
  "EcalEndcapNRemnantClusters",
  "CalorimeterTrackProjections",
  103,
  false,
  {-3.6, -1.7},
  11,
  22
};

const Options ElectronPhotonPEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonPEndCap/ElectronPhotonPEndCap.root",
  "EcalEndcapPExpectedClusters",
  "EcalEndcapPRemnantClusters",
  "CalorimeterTrackProjections",
  102,
  false,
  {1.3, 3.6},
  11,
  22
};

const Options ElectronPhotonBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonBarrel/ElectronPhotonBarrel.root",
  "EcalBarrelExpectedClusters",
  "EcalBarrelRemnantClusters",
  "CalorimeterTrackProjections",
  101,
  false,
  {-1.8, 1.5},
  11,
  22
};

const Options PionNeutronLFOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/PionNeutron_reco.root",
  "Plots/PionNeutronLFHCAL/PionNeutronLFHCAL.root",
  "LFHCALExpectedClusters",
  "LFHCALRemnantClusters",
  "CalorimeterTrackProjections",
  116, 
  false,
  {1.0, 3.5},   
  211,
  2112
};



const Options PionNeutronNEndCapOptions = { 
  "/run/media/epic/Data/Benchmarking/Pairs/PionNeutron_reco.root",
  "Plots/PionNeutronNEndCap/PionNeutronNEndCap.root",
  "HcalEndcapNExpectedClusters",
  "HcalEndcapNRemnantClusters",
  "CalorimeterTrackProjections",
  113, 
  false,
  {-3.5, -1.0}, 
  211,
  2112
};

const Options PionNeutronBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/PionNeutron_reco.root",
  "Plots/PionNeutronBarrel/PionNeutronBarrel.root",
  "HcalBarrelExpectedClusters",
  "HcalBarrelRemnantClusters",
  "CalorimeterTrackProjections",
  111, 
  false,
  {-1.2, 1.2},  
  211,
  2112
};


static double DeltaPhi(double phi1, double phi2) {
  double dphi = phi2 - phi1;
  if (dphi > TMath::Pi()) dphi -= 2.0 * TMath::Pi();
  else if (dphi <= -TMath::Pi()) dphi += 2.0 * TMath::Pi();
  return dphi;
}

// p_track from the projection closest (in eta-phi) to a given position.
// Returns -1 if there is no projection for this system / surface.
static double TrackMomentumAt(const edm4eic::TrackSegmentCollection& segments,
                              const edm4hep::Vector3f& pos, uint64_t system, uint64_t surface) {
  const double eta = edm4hep::utils::eta(pos);
  const double phi = edm4hep::utils::angleAzimuthal(pos);
  double best = numeric_limits<double>::max();
  double p    = -1.0;
  for (edm4eic::TrackSegment segment : segments) {
    for (edm4eic::TrackPoint pr : segment.getPoints()) {
      if (pr.system != system) continue;
      if (system != 101 && pr.surface != surface) continue;
      const double d = std::hypot(eta - edm4hep::utils::eta(pr.position),
                                  DeltaPhi(phi, edm4hep::utils::angleAzimuthal(pr.position)));
      if (d < best) {
        best = d;
        p = std::sqrt(pr.momentum.x * pr.momentum.x + pr.momentum.y * pr.momentum.y +
                      pr.momentum.z * pr.momentum.z);
      }
    }
  }
  return p;
}

void DoubleParticleBenchmark(const Options& opt = PionNeutronNEndCapOptions, uint64_t s_use = 1) {

  gSystem->mkdir(gSystem->DirName(opt.out_file.c_str()), true);
  TFile* Output = new TFile(opt.out_file.c_str(), "RECREATE");

  const int    nDR   = 20;
  const double drMax = 0.5;

  // --- MC truth ---
  TH1D* h_EMC_Charged   = new TH1D("h_EMC_Charged", "E_{MC} charged;E_{MC} [GeV];Counts", 40, 0, 20);
  TH1D* h_EMC_Neutral   = new TH1D("h_EMC_Neutral", "E_{MC} neutral;E_{MC} [GeV];Counts", 40, 0, 20);
  TH1D* h_DeltaR_MC_charged_neutral = new TH1D("h_DeltaR_MC_charged_neutral", "#DeltaR_{MC}(e,#gamma);#DeltaR;Counts", 25, 0, drMax);

  // --- Multiplicities ---
  TH1D* h_NExpected = new TH1D("h_NExpected", "Expected clusters per event;N;Events", 10, -0.5, 9.5);
  TH1D* h_NRemnant  = new TH1D("h_NRemnant",  "Remnant clusters per event;N;Events", 10, -0.5, 9.5);

  // --- Energies ---
  TH1D* h_EExp   = new TH1D("h_EExp", "Expected cluster (charged);E_{exp} [GeV];Events", 40, 0, 20);
  TH1D* h_ERem   = new TH1D("h_ERem", "Remnant clusters (neutral);E_{rem} [GeV];Events", 40, 0, 20);

  // --- Response: reco / MC ---
  TH1D* h_Ratio_Charged = new TH1D("h_Ratio_Charged", "E_{exp} / E_{MC};Ratio;Events", 50, 0, 2);
  TH1D* h_Ratio_Neutral = new TH1D("h_Ratio_Neutral", "E_{rem} / E_{MC};Ratio;Events", 50, 0, 2);
  TH1D* h_Ratio_Tot  = new TH1D("h_Ratio_Tot", "(E_{exp}+E_{rem}) / (E_{MC-pair});Ratio;Events", 50, 0, 2);

  // --- Comparison to the track ---
  TH1D* h_ExpOverPTrack  = new TH1D("h_ExpOverPTrack", "E_{exp} / p_{track};Ratio;Events", 100, 0, 2);
  TH1D* h_DeltaR_MC_Charged = new TH1D("h_DeltaR_MC_Charged", "Cluster distance to MC charged;#DeltaR;Events", 100, 0, 2);
  TH1D* h_DeltaR_MC_Neutral = new TH1D("h_DeltaR_MC_Neutral", "Cluster distance to MC neutral;#DeltaR;Events", 100, 0, 2);


  // --- 2D: reco vs MC ---
  TH2D* h2_EExp_vs_EMC = new TH2D("h2_EExp_vs_EMC", "Charged: E_{exp} vs E_{MC};E_{MC} [GeV];E_{exp} [GeV]", 40, 0, 20, 40, 0, 20);
  TH2D* h2_ERem_vs_EMC = new TH2D("h2_ERem_vs_EMC", "Neutral: E_{rem} vs E_{MC};E_{MC} [GeV];E_{rem} [GeV]", 40, 0, 20, 40, 0, 20);

  // --- Dependence on the charged-neutral distance (MC level) ---
  TProfile* p_RatioCharged_vs_dR = new TProfile("p_RatioCharged_vs_dR", "<E_{exp}/E_{MC}> vs #DeltaR_{MC};#DeltaR_{MC};<E_{exp}/E_{MC}>", nDR, 0, drMax);
  TProfile* p_RatioNeutral_vs_dR = new TProfile("p_RatioNeutral_vs_dR", "<E_{rem}/E_{MC}> vs #DeltaR_{MC};#DeltaR_{MC};<E_{rem}/E_{MC}>", nDR, 0, drMax);
  TProfile* p_RatioCharged_vs_Eta = new TProfile("p_RatioCharged_vs_Eta", "<E_{exp}/E_{MC}> ;#eta;<E_{exp}/E_{MC}>", 30, opt.etaRange[0], opt.etaRange[1]);
  TProfile* p_RatioCharged_vs_E = new TProfile("p_RatioCharged_vs_E", "<E_{exp}/E_{MC}> ;E_{MC} [GeV];<E_{exp}/E_{MC}>", 40, 0, 20);
  TProfile* p_RatioNeutral_vs_Eta = new TProfile("p_RatioNeutral_vs_Eta", "<E_{rem}/E_{MC}> ;#eta;<E_{rem}/E_{MC}>", 30, opt.etaRange[0], opt.etaRange[1]);
  TProfile* p_RatioNeutral_vs_E = new TProfile("p_RatioNeutral_vs_E", "<E_{rem}/E_{MC}> ;E_{MC} [GeV];<E_{rem}/E_{MC}>", 40, 0, 20);

  TProfile* p_RatioTot_vs_dR  = new TProfile("p_RatioTot_vs_dR", "<E_{tot}^{reco}/E_{tot}^{MC}> vs #DeltaR_{MC};#DeltaR_{MC};ratio", nDR, 0, drMax);
  TProfile* p_ExpOverP_vs_dR  = new TProfile("p_ExpOverP_vs_dR", "<E_{exp}/p_{track}> vs #DeltaR_{MC};#DeltaR_{MC};<E_{exp}/p_{track}>", nDR, 0, drMax);
  TH2D* h2_RatioNeutral_vs_dR = new TH2D("h2_RatioNeutral_vs_dR", "E_{rem}/E_{MC} vs #DeltaR_{MC};#DeltaR_{MC};E_{rem}/E_{MC}", nDR, 0, drMax, 100, 0, 2);
  TH2D* h2_RatioCharged_vs_dR = new TH2D("h2_RatioCharged_vs_dR", "E_{exp}/E_{MC} vs #DeltaR_{MC};#DeltaR_{MC};E_{exp}/E_{MC}", nDR, 0, drMax, 100, 0, 2);

  // --- Reconstruction efficiency vs distance / energy ---
  // "found" = at least one Expected (charged) / Remnant (neutral) cluster
  TEfficiency* eff_Charged_vs_dR = new TEfficiency("eff_Charged_vs_dR", "Charged (expected cluster) found;#DeltaR_{MC};Efficiency", nDR, 0, drMax);
  TEfficiency* eff_Neutral_vs_dR = new TEfficiency("eff_Neutral_vs_dR", "Neutral (remnant cluster) found;#DeltaR_{MC};Efficiency", nDR, 0, drMax);
  TEfficiency* eff_Charged_vs_E  = new TEfficiency("eff_Charged_vs_E", "Charged (expected cluster) found;E_{MC} [GeV];Efficiency", 40, 0, 20);
  TEfficiency* eff_Neutral_vs_E  = new TEfficiency("eff_Neutral_vs_E", "Neutral (remnant cluster) found;E_{MC} [GeV];Efficiency", 40, 0, 20);
  TEfficiency* eff_Charged_vs_Eta  = new TEfficiency("eff_Charged_vs_Eta", "Charged (expected cluster) found;#eta;Efficiency", 30, opt.etaRange[0], opt.etaRange[1]);
  TEfficiency* eff_Neutral_vs_Eta  = new TEfficiency("eff_Neutral_vs_Eta", "Neutral (remnant cluster) found;#eta;Efficiency", 30, opt.etaRange[0], opt.etaRange[1]);

  cout << "\n  Start macro: Expected/Remnant cluster analysis" << endl;

  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile(opt.in_file);

  uint64_t nFrames = reader.getEntries(podio::Category::Event);
  if (opt.debug) nFrames = std::min<uint64_t>(nFrames, 20);
  cout << "  Processing " << nFrames << " frames..." << endl;

  uint64_t nUsed = 0, nChargedFound = 0, nNeutralFound = 0;

  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {
    auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));

    auto& expected     = frame.get<edm4eic::ClusterCollection>(opt.expectedClusters);
    auto& remnant      = frame.get<edm4eic::ClusterCollection>(opt.remnantClusters);
    auto& segments     = frame.get<edm4eic::TrackSegmentCollection>(opt.projections);
    auto& mc_particles = frame.get<edm4hep::MCParticleCollection>("MCParticles");

    // 1. Charged and neutral MC truth
    TLorentzVector mcCharged, mcNeutral;
    bool foundCharged = false, foundNeutral = false;

    bool first=true;

    for (edm4hep::MCParticle mc_part : mc_particles) {
      if (mc_part.getGeneratorStatus() != 1) continue;
      const int pdg = mc_part.getPDG();
      auto mom = mc_part.getMomentum();
      const double eta_mc = edm4hep::utils::eta(mom);
      if (eta_mc < opt.etaRange[0] || eta_mc > opt.etaRange[1]) continue;
      if(first) {
        h_NExpected->Fill(expected.size());
        h_NRemnant->Fill(remnant.size());
        first=false;
      }
      

      if (!foundCharged && std::abs(pdg) == std::abs(opt.pdgExpected)) {
        mcCharged.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
        foundCharged = true;
      } else if (!foundNeutral && std::abs(pdg) == std::abs(opt.pdgRemnant)) {
        mcNeutral.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
        foundNeutral = true;
      }
    }
    if (!foundCharged || !foundNeutral) {
      if (opt.debug) cout << "Frame " << iFrame << ": no charged-neutral MC pair. Skipping." << endl;
      continue;
    }
    ++nUsed;

    const double eMCCharged = mcCharged.E();
    const double etaMCCharged = mcCharged.Eta();

    const double eMCNeutral = mcNeutral.E();
    const double etaMCNeutral = mcNeutral.Eta();

    const double dRMC = std::hypot(etaMCCharged - etaMCNeutral, DeltaPhi(mcCharged.Phi(), mcNeutral.Phi()));

    h_EMC_Charged->Fill(eMCCharged);
    h_EMC_Neutral->Fill(eMCNeutral);
    h_DeltaR_MC_charged_neutral->Fill(dRMC);

    // 2. Charged = Expected cluster (the most energetic if there are several)
    bool   hasExp = false;
    double eExp = 0.0;
    double DR_Exp_track=0.0;
    edm4hep::Vector3f posExp{0, 0, 0};
    for (edm4eic::Cluster c : expected) {
      posExp = c.getPosition();
      DR_Exp_track=std::hypot(etaMCCharged - edm4hep::utils::eta(posExp), DeltaPhi(mcCharged.Phi(), edm4hep::utils::angleAzimuthal(posExp)));
      h_DeltaR_MC_Charged->Fill(DR_Exp_track);
      if(DR_Exp_track<0.6){
        if (!hasExp || c.getEnergy() > eExp) {
          eExp = c.getEnergy();
          hasExp = true;
        }
      }
    }

    // 3. Neutral = Remnant cluster leading
    bool   hasRem = false;
    double eRem = 0.0;
    double DR_Rem_track=0.0;
    edm4hep::Vector3f posRem{0, 0, 0};
    for (edm4eic::Cluster c : remnant) {
      posRem = c.getPosition();
      DR_Rem_track=std::hypot(etaMCNeutral - edm4hep::utils::eta(posRem), DeltaPhi(mcNeutral.Phi(), edm4hep::utils::angleAzimuthal(posRem)));
      h_DeltaR_MC_Neutral->Fill(DR_Rem_track);
      if(DR_Rem_track<0.6){
        if (!hasRem || c.getEnergy() > eRem) {
          eRem = c.getEnergy();
          hasRem = true;
        }
      }
    }
    hasExp = hasExp && eExp > 0;
    hasRem = hasRem && eRem > 0;

    // 4. Efficiencies (denominator = all charged-neutral MC events)
    eff_Charged_vs_dR->Fill(hasExp, dRMC);
    eff_Neutral_vs_dR->Fill(hasRem, dRMC);
    eff_Charged_vs_E->Fill(hasExp, eMCCharged);
    eff_Neutral_vs_E->Fill(hasRem, eMCNeutral);
    eff_Charged_vs_Eta->Fill(hasExp, etaMCCharged);
    eff_Neutral_vs_Eta->Fill(hasRem, etaMCNeutral);
    if (hasExp) ++nChargedFound;
    if (hasRem) ++nNeutralFound;

    // 5. Charged energy response + comparison to the track
    double pTrack = -1.0;
    if (hasExp) {
      h_EExp->Fill(eExp);
      h_Ratio_Charged->Fill(eExp / eMCCharged);
      h2_EExp_vs_EMC->Fill(eMCCharged, eExp);
      p_RatioCharged_vs_dR->Fill(dRMC, eExp / eMCCharged);
      p_RatioCharged_vs_Eta->Fill(etaMCCharged, eExp / eMCCharged);
      p_RatioCharged_vs_E->Fill(eMCCharged, eExp/ eMCCharged);
      h2_RatioCharged_vs_dR->Fill(dRMC, eExp / eMCCharged);

      pTrack = TrackMomentumAt(segments, posExp, opt.system, s_use);
      if (pTrack > 0) {
        h_ExpOverPTrack->Fill(eExp / pTrack);
        p_ExpOverP_vs_dR->Fill(dRMC, eExp / pTrack);
      }
    }

    // 6. Neutral energy response
    if (hasRem) {
      h_ERem->Fill(eRem);
      h_Ratio_Neutral->Fill(eRem / eMCNeutral);
      h2_ERem_vs_EMC->Fill(eMCNeutral, eRem);
      p_RatioNeutral_vs_dR->Fill(dRMC, eRem / eMCNeutral);
      p_RatioNeutral_vs_Eta->Fill(etaMCNeutral, eRem / eMCNeutral);
      p_RatioNeutral_vs_E->Fill(eMCNeutral, eRem / eMCNeutral);
      h2_RatioNeutral_vs_dR->Fill(dRMC, eRem / eMCNeutral);
    }

    // 7. Remnant <-> Expected distance in eta-phi (reco level)
    if (hasExp && hasRem) {
      const double ratioTot = (eExp + eRem) / (eMCCharged + eMCNeutral);
      h_Ratio_Tot->Fill(ratioTot);
      p_RatioTot_vs_dR->Fill(dRMC, ratioTot);
    }

    if (opt.debug) {
      cout << "Frame " << iFrame << ": dR_MC=" << dRMC
           << " | E_MC(charged,neutral)=" << eMCCharged << "," << eMCNeutral
           << " | E_exp=" << eExp << " | E_rem()=" << eRem
           << " | p_track=" << pTrack << endl;
    }
  }  // end event loop

  cout << "  Events with charged-neutral MC pair: " << nUsed << endl;
  cout << "  Expected (charged) found:   " << nChargedFound << "/" << nUsed << endl;
  cout << "  Remnant (neutral) found:    " << nNeutralFound << "/" << nUsed << endl;

  Output->Write();
  Output->Close();

  cout << "  Analysis complete! Histograms saved to: " << opt.out_file << endl;
}