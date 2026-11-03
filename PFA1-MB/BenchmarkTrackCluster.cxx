// ----------------------------------------------------------------------------
//
// Quick ROOT macro to check matching calorimeter clusters with track projections.
//
// Usage:
//   root -b -q BenchmarkTrackCluster.cxx++
// ----------------------------------------------------------------------------

#define MATCHPROJECTIONSANDCLUSTERS_CXX

// c++ utilities
#include <cmath>
#include <string>
#include <limits>
#include <optional>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>
// root libraries
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
// podio libraries
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/CollectionBase.h>
// edm4eic types
#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/TrackPoint.h>
#include <edm4eic/TrackSegment.h>
#include <edm4eic/TrackSegmentCollection.h>
// edm4hep types
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticleCollection.h>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::pair;
using std::optional;
using std::numeric_limits;

// user options ---------------------------------------------------------------

struct Options {
  string in_file;                   // input file
  string out_file;                  // output file
  string clusters;                  // name of collection of clusters to match to
  string projections;               // name of collection of track projections
  uint64_t system;                  // system id of calo to match to
  bool debug;                       // enables info messages
  double MaxREtaPhi;             // sets the boundary of good Track-Cluster distance
  bool find_DR;
  std::array<double, 2> etaRange;   // [eta_min, eta_max] for histogram drawing
  std::array<double, 2> ptRange;    // [pt_min, pt_max] for histogram drawing
    // ECAL info (only for pion files, leave empty to disable)
  string ecal_clusters;             // name of ECAL cluster collection
  uint64_t ecal_system;             // system id of ECAL
  double ecal_MaxREtaPhi;           // track-ECAL cluster matching cut
};

const Options ElectronPEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Electron/Electron.PEndCap.root",
  "Plots/ElectronPEndCap.root",
  "EcalEndcapPClusters",
  "CalorimeterTrackProjections",
  102, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04,
  false,
  {1.3, 3.6},  
  {0.0, 10.0}  
};

const Options ElectronNEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Electron/Electron_2.NEndCap.root",
  "Plots/ElectronNEndCap.root",
  "EcalEndcapNClusters",
  "CalorimeterTrackProjections",
  103, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04,
  false,
  {-3.6, -1.8}, 
  {0.0, 7.0}   
};

const Options ElectronBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Electron/Electron_2.Barrel.root",
  "Plots/ElectronBarrel.root",
  "EcalBarrelClusters",
  "CalorimeterTrackProjections",
  101, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04,
  false,
  {-1.8, 1.5},  
  {0.0, 20.0}   
};

const Options PionLFOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Pion/Pion_2.LFHCAL.root",
  "Plots/PionLFHCAL.root",
  "LFHCALClusters",
  "CalorimeterTrackProjections",
  116, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.1,
  false,
  {1.0, 3.5},   
  {0.0, 13.0},   
  "EcalEndcapPClusters", 102, 0.04
};

const Options PionPEndCapOptions = { 
  "/run/media/epic/Data/Benchmarking/Single/Pion/Pion_2.Insert.root",
  "Plots/PionPEndCap.root",
  "HcalEndcapPInsertClusters",
  "CalorimeterTrackProjections",
  116, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.06,
  false,
  {3.1, 3.8},   
  {0.0, 2.0}   
};

const Options PionNEndCapOptions = { 
  "/run/media/epic/Data/Benchmarking/Single/Pion/Pion_2.NEndCap.root",
  "Plots/PionNEndCap.root",
  "HcalEndcapNClusters",
  "CalorimeterTrackProjections",
  113, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.2,
  false,
  {-3.5, -1.0}, 
  {0.0, 13.0},
  "EcalEndcapNClusters", 103, 0.04  
};

const Options PionBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Pion/Pion_2.Barrel.root",
  "Plots/PionBarrel.root",
  "HcalBarrelClusters",
  "CalorimeterTrackProjections",
  111, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.06,
  false,
  {-1.2, 1.2},  
  {0.0, 20.0},
  "EcalBarrelClusters", 101, 0.04  
};

double DeltaPhi(double phi1, double phi2) {
    double dphi = phi2 - phi1;
    if (dphi > TMath::Pi()) {
        dphi -= 2.0 * TMath::Pi();
    } else if (dphi <= -TMath::Pi()) {
        dphi += 2.0 * TMath::Pi();
    }
    return dphi;
}

// macro body -----------------------------------------------------------------
void BenchmarkTrackCluster(const Options& opt = PionLFOptions, uint64_t s_use = 1) {

    TFile *Output = new TFile(opt.out_file.c_str(), "RECREATE");

    // Read the provided pt and eta ranges from the opt structure
    const double ptMin  = opt.ptRange[0];
    const double ptMax  = opt.ptRange[1];
    const double etaMin = opt.etaRange[0];
    const double etaMax = opt.etaRange[1];

    // MC Particles histograms 
    TH1D* MCMomentum      = new TH1D("MCMomentum", "MCMomentum;p [GeV];Counts", 40, 0, 20);
    TH1D* MCEta           = new TH1D("MCEta", "MCEta;#eta;Counts", 40, etaMin, etaMax);
    TH1D* MCPt            = new TH1D("MCPt", "MCPt;p_{T} [GeV];Counts", 40, ptMin, ptMax);
    TH1D* MCPhi           = new TH1D("MCPhi", "MCPhi;#phi [rad];Counts", 40, -3.15, 3.15);

    TH1D* MCMomentumHited = new TH1D("MCMomentumHited", "MCMomentumHited;p [GeV];Counts", 40, 0, 20);
    TH1D* MCEtaHited      = new TH1D("MCEtaHited", "MCEtaHited;#eta;Counts", 40, etaMin, etaMax);
    TH1D* MCPtHited       = new TH1D("MCPtHited", "MCPtHited;p_{T} [GeV];Counts", 40, ptMin, ptMax);
    TH1D* MCPhiHited      = new TH1D("MCPhiHited", "MCPhiHited;#phi [rad];Counts", 40, -3.15, 3.15);

    // Track-Cluster matching histograms 
    TH1D* REtaPhiHist     = new TH1D("REtaPhiHist", "REtaPhiHist;r;Events", 50, 0, opt.MaxREtaPhi*5);
    TH2D* DeltaEtaPhiHist = new TH2D("DeltaEtaPhiHist", "DeltaEtaPhiHist;#Delta#eta;#Delta#phi", 100, -opt.MaxREtaPhi*2, opt.MaxREtaPhi*2, 100, -opt.MaxREtaPhi*2, opt.MaxREtaPhi*2);

    // ECluster/ETrack histograms 
    TH1D* EClusteroverETrack = new TH1D("EClusteroverETrack", "E_{Cl}/E_{Tr};Ratio E_{Cl}/E_{Tr};Events", 50, 0, 2);
    TH1D* EMCoverECluster    = new TH1D("EMCoverECluster", "E_{MC}/E_{Cl};Ratio E_{MC}/E_{Cl};Events", 50, 0, 2);
    TH1D* EMCoverETrack      = new TH1D("EMCoverETrack", "E_{MC}/E_{Tr};Ratio E_{MC}/E_{Tr};Events", 50, 0, 2);

    TH2D* EnergyRatiovsTrackPt  = new TH2D("EnergyRatiovsTrackPt", "E_{Cl}/E_{Tr} vs TrackPt;p_{T} [GeV];Ratio E_{Cl}/E_{Tr}", 40, ptMin, ptMax, 50, 0, 2);
    TH2D* EnergyRatiovsTrackEta = new TH2D("EnergyRatiovsTrackEta", "EnergyRatiovsTrackEta;#eta;Ratio E_{Cl}/E_{Tr}", 40, etaMin, etaMax, 50, 0, 2);

    // Multiplicity of tracks and clusters
    TH1D* NClustersPerEvent      = new TH1D("NClustersPerEvent", "Number of clusters per event;N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* NClustersMatchedToTrack = new TH1D("NClustersMatchedToTrack", "Number of matched clusters per track;N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* NTracksMatchedToCluster = new TH1D("NTracksMatchedToCluster", "Number of matched tracks per cluster;N_{tracks};Events", 10, -0.5, 9.5);

    // ECluster-f*ETrack
    TH1D* EClustMinusFsubEtrk    = new TH1D("EClustMinusFsubEtrk", "E_{clust} - f_{sub}E_{trk};E_{clust} - f_{sub}E_{trk} [GeV];Events", 40, -10, 10);
 
    // ECAL+HCAL to track (pion files only): exactly 4 histograms
    const bool doEcal = !opt.ecal_clusters.empty();
    TH1D* EEcalHcalOverETrack  = new TH1D("EEcalHcalOverETrack", "(E_{ECal}+E_{HCal})/E_{Tr};Ratio (E_{ECal}+E_{HCal})/E_{Tr};Events", 50, 0, 2);
    // announce start of macro
    cout << "\n  Beginning cluster-track projection matching macro!" << endl;

    // open file w/ frame reader
    podio::ROOTReader reader = podio::ROOTReader();
    reader.openFile(opt.in_file);
    cout << "    Opened ROOT-based frame reader." << endl;
        
    // get no. of frames and announce
    uint64_t nFrames = reader.getEntries(podio::Category::Event);
    if (opt.debug) nFrames = 10;
    cout << "    Starting frame loop: " << nFrames << " frames to process." << endl;

    // iterate through frames (i.e. events in this case)
    uint64_t nClustTotal   = 0;
    uint64_t nClustMatched = 0;

    vector<double> mcEtaAll;
    vector<double> mcEtaHit;

    vector<double> vClustE;
    vector<double> vPTrack;
    vector<double> vDistMatch;

    for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {
        // announce progress
        if (opt.debug) cout << "      Processing frame " << iFrame + 1 << "/" << nFrames << "..." << endl;

        // grab frame
        auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));

        // grab collections
        auto& clusters     = frame.get<edm4eic::ClusterCollection>(opt.clusters);
        auto& segments     = frame.get<edm4eic::TrackSegmentCollection>(opt.projections);
        auto& mc_particles = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        const edm4eic::ClusterCollection* ecalClusters = nullptr;
        if (doEcal) ecalClusters = &frame.get<edm4eic::ClusterCollection>(opt.ecal_clusters);

        NClustersPerEvent->Fill(clusters.size());
        if (opt.debug) cout << "      Number of MCParticles in frame: " << mc_particles.size() << endl;
        
        // loop over MCParticle
        TLorentzVector Particle;
        for (edm4hep::MCParticle mc_part : mc_particles) {
            if (mc_part.getGeneratorStatus() != 1) continue;

            int pdg       = mc_part.getPDG();              
            auto momentum = mc_part.getMomentum();   
            
            Particle.SetPxPyPzE(momentum.x, momentum.y, momentum.z, mc_part.getEnergy());
            
            mcEtaAll.push_back(Particle.Eta());

            MCMomentum->Fill(Particle.P());
            MCEta->Fill(Particle.Eta());
            MCPt->Fill(Particle.Pt());
            MCPhi->Fill(Particle.Phi());

            if (opt.debug) {
                cout << "      [MC Particle] PDG: " << pdg 
                     << " | E: " << Particle.E() << " GeV"
                     << " | Eta: " << Particle.Eta() 
                     << " | Phi: " << Particle.Phi() << endl;
            }
        }

        map<int32_t, int> trackMatchCount;
        bool particleCountedAsHit = false;

        // loop over clusters
        for (edm4eic::Cluster cluster : clusters) {

            int tracksInsideCut = 0;

            // grab eta/phi of cluster
            const double etaClust = edm4hep::utils::eta(cluster.getPosition());
            const double phiClust = edm4hep::utils::angleAzimuthal(cluster.getPosition());

            // match based on eta/phi distance
            double distMatch = numeric_limits<double>::max(); 
            double dEta      = numeric_limits<double>::max(); 
            double dPhi      = numeric_limits<double>::max(); 

            // loop over projections to find matching one
            optional<edm4eic::TrackPoint> match;
            optional<edm4eic::TrackSegment> matchedSegment;
            int SegID = 0, BestSegID = -1;

            for (edm4eic::TrackSegment segment : segments) {
                for (edm4eic::TrackPoint projection : segment.getPoints()) {
                    if (opt.debug) cout << " Projection-system: " << projection.system << "|  opt.system " << opt.system << "| surface: " << projection.surface << endl;
                    
                    const bool isInSystem = (projection.system == opt.system);
                    const bool isAtFace   = (projection.surface == s_use);

                    if (opt.system == 101) {
                        if (!isInSystem) continue;
                    } else {
                        if (!isInSystem || !isAtFace) continue;
                    }

                    // grab eta/phi of projection
                    const double etaProject = edm4hep::utils::eta(projection.position);
                    const double phiProject = edm4hep::utils::angleAzimuthal(projection.position);
                    if (opt.debug) cout << "Eta Track: " << etaProject << endl;

                    // get distance to projection
                    const double distProject = hypot(
                        etaClust - etaProject,
                        DeltaPhi(phiClust, phiProject)
                    );

                    // if smallest distance found, update variables accordingly
                    if (distProject < distMatch) {
                        distMatch      = distProject;
                        dEta           = etaClust - etaProject;
                        dPhi           = DeltaPhi(phiClust, phiProject);
                        match          = projection;
                        matchedSegment = segment;
                        BestSegID      = SegID;
                        tracksInsideCut++;
                    }
                }  // end point loop
                SegID++;
            }  // end segment loop 

            NTracksMatchedToCluster->Fill(tracksInsideCut);

            DeltaEtaPhiHist->Fill(dEta, dPhi);
            REtaPhiHist->Fill(distMatch);
            ++nClustTotal;

            // do analysis if match was found...
            if (match.has_value() && distMatch < opt.MaxREtaPhi) {
                edm4eic::TrackPoint matchedProject = match.value();
                if (opt.debug) cout << "        Found match! match = " << matchedProject << endl;
                
                auto TrackMomentum = matchedProject.momentum;
                double pTrack = sqrt(TrackMomentum.x * TrackMomentum.x + 
                                     TrackMomentum.y * TrackMomentum.y + 
                                     TrackMomentum.z * TrackMomentum.z);
                
                // High energy limit pTrack=eTrack
                const double EnergyRatio = cluster.getEnergy() / pTrack;
                EClusteroverETrack->Fill(EnergyRatio);
                EMCoverECluster->Fill(cluster.getEnergy() / Particle.E());
                EMCoverETrack->Fill(pTrack / Particle.E());

                const double pT_Track = hypot(TrackMomentum.x, TrackMomentum.y);
                EnergyRatiovsTrackPt->Fill(pT_Track, EnergyRatio);
                
                const double etaProject = edm4hep::utils::eta(matchedProject.position);
                EnergyRatiovsTrackEta->Fill(etaProject, EnergyRatio);

                if (doEcal) {
                    // find the ECAL projection of the SAME track
                    double distEcal = numeric_limits<double>::max();
                    double eEcal    = 0.0;
                    for (edm4eic::TrackPoint pt : matchedSegment->getPoints()) {
                        if (pt.system != opt.ecal_system) continue;
                        if (opt.ecal_system != 101 && pt.surface != s_use) continue;

                        const double etaPE = edm4hep::utils::eta(pt.position);
                        const double phiPE = edm4hep::utils::angleAzimuthal(pt.position);

                        // closest ECAL cluster to this projection
                        for (edm4eic::Cluster ecl : *ecalClusters) {
                            const double d = hypot(
                                edm4hep::utils::eta(ecl.getPosition()) - etaPE,
                                DeltaPhi(edm4hep::utils::angleAzimuthal(ecl.getPosition()), phiPE)
                            );
                            if (d < distEcal && d < opt.ecal_MaxREtaPhi) {
                                distEcal = d;
                                eEcal    = ecl.getEnergy();
                            }
                        }
                    }

                    const double rHcal = cluster.getEnergy() / pTrack;
                    const double rEcal = eEcal / pTrack;
                    EEcalHcalOverETrack->Fill(rEcal + rHcal);
                }

                ++nClustMatched;
                
                trackMatchCount[BestSegID]++;

                if (!particleCountedAsHit) {
                    mcEtaHit.push_back(Particle.Eta());
                    particleCountedAsHit = true;
                    MCMomentumHited->Fill(Particle.P());
                    MCEtaHited->Fill(Particle.Eta());
                    MCPtHited->Fill(Particle.Pt());
                    MCPhiHited->Fill(Particle.Phi());
                }

                // stash values needed to auto-compute f_sub and thresh below
                vClustE.push_back(cluster.getEnergy());
                vPTrack.push_back(pTrack);
                vDistMatch.push_back(distMatch);
            }
        }  // end cluster loop

        for (const auto& [trackID, count] : trackMatchCount) {
            if (trackID != -1) NClustersMatchedToTrack->Fill(count);
        }
    }  
    // end frame loop

    if (!mcEtaHit.empty()) {
        double minEtaFound = *std::min_element(mcEtaHit.begin(), mcEtaHit.end());
        double maxEtaFound = *std::max_element(mcEtaHit.begin(), mcEtaHit.end());

        int nInRange = std::count_if(mcEtaAll.begin(), mcEtaAll.end(),
                                     [minEtaFound, maxEtaFound](double eta) {
                                         return eta >= minEtaFound && eta <= maxEtaFound;
                                     });

        cout << "  Eta Range found: [" << minEtaFound << ", " << maxEtaFound << "]" << endl;
        cout << "  Number of MCTruth in this range: " << nInRange << endl;
        
        cout << "  Particles matched " << MCEtaHited->GetEntries() << "/" << nInRange << " particles in range" << endl;
        if (MCEtaHited->GetEntries() > 0) {
            cout << "  Average number of cluster matched " << (double)nClustMatched / MCEtaHited->GetEntries() << " to particle" << endl;
        }
    }

    cout << "    Finished frame loop: " << nClustMatched << "/" << nClustTotal << " clusters matched." << endl;

    const double f_sub = EClusteroverETrack->GetMean();
    cout << "    Auto-computed f_sub = " << f_sub << endl;

    for (size_t i = 0; i < vClustE.size(); ++i) {
        EClustMinusFsubEtrk->Fill(vClustE[i] - f_sub * vPTrack[i]);
    }

    if (opt.find_DR) {
        std::vector<double> drCuts = {0.02, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
        const double E_MIN_THRESHOLD = 0.1; 

        std::cout << "\n  --- dR cut scan ---\n";
        std::cout << std::left
                  << std::setw(8)  << "dR"
                  << std::setw(10) << "Nmatch"
                  << std::setw(12) << "purity"
                  << std::setw(10) << "Ngood"
                  << std::setw(10) << "Nbad"
                  << std::setw(12) << "dGood"
                  << std::setw(12) << "dBad"
                  << std::setw(12) << "good:bad"
                  << std::endl;

        int prevGood = 0;
        int prevBad  = 0;

        for (double drCut : drCuts) {
            double nMatched   = 0;
            double nLowEnergy = 0;

            for (size_t i = 0; i < vDistMatch.size(); ++i) {
                if (vDistMatch[i] < drCut) {
                    nMatched++;
                    if (vClustE[i] / vPTrack[i] < E_MIN_THRESHOLD) nLowEnergy++;
                }
            }

            const int nGood = nMatched - nLowEnergy;
            const int nBad  = nLowEnergy;
            const double purity = (nMatched > 0) ? double(nGood) / nMatched : 0.0;

            const int dGood = nGood - prevGood;
            const int dBad  = nBad  - prevBad;
            const double marginalRatio = (dBad > 0) ? double(dGood) / dBad : std::numeric_limits<double>::infinity();

            std::cout << std::left
                      << std::setw(8)  << drCut
                      << std::setw(10) << nMatched
                      << std::setw(12) << purity
                      << std::setw(10) << nGood
                      << std::setw(10) << nBad
                      << std::setw(12) << dGood
                      << std::setw(12) << dBad
                      << std::setw(12) << marginalRatio
                      << std::endl;

            prevGood = nGood;
            prevBad  = nBad;
        }
        std::cout << "  --- end of dR cut scan ---\n" << std::endl;
    }

    TProfile* EnergyRatioProfilevsTrackPt = EnergyRatiovsTrackPt->ProfileX();
    EnergyRatioProfilevsTrackPt->SetName("EnergyRatioProfilevsTrackPt");
    EnergyRatioProfilevsTrackPt->SetTitle("Mean E_{Cl}/E_{Tr} vs Track p_{T};p_{T} [GeV];<E_{Cl}/E_{Tr}>");

    TProfile* EnergyRatioProfilevsTrackEta = EnergyRatiovsTrackEta->ProfileX();
    EnergyRatioProfilevsTrackEta->SetName("EnergyRatioProfilevsTrackEta");
    EnergyRatioProfilevsTrackEta->SetTitle("Mean E_{Cl}/E_{Tr} vs Track Eta;#eta;<E_{Cl}/E_{Tr}>");

    MCMomentum->Write();
    MCEta->Write();
    MCPt->Write();
    MCPhi->Write();
    MCMomentumHited->Write();
    MCEtaHited->Write();
    MCPtHited->Write();
    MCPhiHited->Write();

    DeltaEtaPhiHist->Write();
    REtaPhiHist->Write();

    EClusteroverETrack->Write();
    EMCoverECluster->Write();
    EMCoverETrack->Write();

    EnergyRatiovsTrackPt->Write();
    EnergyRatioProfilevsTrackPt->Write();

    EnergyRatiovsTrackEta->Write();
    EnergyRatioProfilevsTrackEta->Write();

    NClustersPerEvent->Write();
    NTracksMatchedToCluster->Write();
    NClustersMatchedToTrack->Write();

    EClustMinusFsubEtrk->Write();

    if (doEcal) {
        EEcalHcalOverETrack->Write();
    }

    Output->Close();
    
    // announce end and exit
    cout << "  End of macro!\n" << endl;
    return;
}

// end ------------------------------------------------------------------------