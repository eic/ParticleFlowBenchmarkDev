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



// user options ---------------------------------------------------------------

struct Options {
  std::string in_file;      // input file
  std::string out_file;     // output file
  std::string clusters;     // name of collection of clusters to match to
  std::string projections;  // name of collection of track projections
  uint64_t    system;       // system id of calo to match to
  bool debug;               // enables info messages
  double MaxREtaPhi;        // sets the boundry of good Track-Cluster distance
};
const Options ElectronOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Electron/*.root",
  "Plots/Electron.root",
  "EcalEndcapPClusters",
  "CalorimeterTrackProjections",
  102, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04
};

const Options PionOptions = {
  "/run/media/epic/Data/Benchmarking/Single/Pion/*.root",
  "Plots/Pion.root",
  "LFHCALClusters",
  "CalorimeterTrackProjections",
   116, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.2
};
double DeltaPhi( double phi1, double phi2) {
    double dphi = phi2-phi1;
    if (dphi > TMath::Pi()) {
        dphi -= 2.0 * TMath::Pi();
    } else if (dphi <= -TMath::Pi()) {
        dphi += 2.0 * TMath::Pi();
    }
    return dphi;
}


// macro body -----------------------------------------------------------------

void BenchmarkTrackCluster(const Options& opt = PionOptions) {


    TFile *Output = new TFile(opt.out_file.c_str(), "RECREATE");
    // MC Particles histograms 
    TH1D* MCMomentum= new TH1D("MCMomentum","MCMomentum",100,-0.5,20.5);
    TH1D* MCEta= new TH1D("MCEta","MCEta",100,0,3.8);
    TH1D* MCPt= new TH1D("MCPt","MCPt",100,-0.5,20.5);
    TH1D* MCPhi= new TH1D("MCPhi","MCPhi",100,-3.15,3.15);

    TH1D* MCMomentumHited= new TH1D("MCMomentumHited","MCMomentumHited",100,-0.5,20.5);
    TH1D* MCEtaHited= new TH1D("MCEtaHited","MCEtaHited",100,0,3.8);
    TH1D* MCPtHited= new TH1D("MCPtHited","MCPtHited",100,-0.5,20.5);
    TH1D* MCPhiHited= new TH1D("MCPhiHited","MCPhiHited",100,-3.15,3.15);

    // Track-Cluster matching histograms 
    TH1D* REtaPhiHist= new TH1D("REtaPhiHist","REtaPhiHist;r;Events",1000,0,0.5);
    //TH2D* DeltaEtaPhiHist= new TH2D("DeltaEtaPhiHist","DeltaEtaPhiHist;#eta,#phi",1000,-0.05,0.05,1000,-0.05,0.05);
    TH2D* DeltaEtaPhiHist= new TH2D("DeltaEtaPhiHist","DeltaEtaPhiHist;#eta,#phi",1000,-0.3,0.3,1000,-0.3,0.3);

    // ECluster/ETrack histograms 
    TH1D* EClusteroverETrack= new TH1D("EClusteroverETrack","E_{Cl}/E_{Tr};Ratio E_{Cl}/E_{Tr};Events",100,0,5);
    TH1D* EMCoverECluster= new TH1D("EMCoverECluster","E_{MC}/E_{Cl};Ratio E_{MC}/E_{Cl};Events",100,0,5);
    TH1D* EMCoverETrack= new TH1D("EMCoverETrack","E_{MC}/E_{Tr};Ratio E_{MC}/E_{Tr};Events",100,0,5);

    TH2D* EnergyRatiovsTrackPt= new TH2D("EnergyRatiovsTrackPt","E_{Cl}/E_{Tr} vs TrackPt;pT [GeV];Ratio E_{Cl}/E_{Tr}",45,0,9,100,0,5);
    TH2D* EnergyRatiovsTrackEta= new TH2D("EnergyRatiovsTrackEta","EnergyRatiovsTrackEta;#eta;Ratio E_{Cl}/E_{Tr}",50,1.25,3.5,100,0,5);

    // Multiplicity of tracks and clusters
    TH1D* NClustersPerEvent = new TH1D("NClustersPerEvent", "Number of clusters per event;N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* NClustersMatchedToTrack = new TH1D("NClustersMatchedToTrack", "Number of matched clusters per track;N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* NTracksMatchedToCluster = new TH1D("NTracksMatchedToCluster", "Number of matched tracks per cluster;N_{tracks};Events", 10, -0.5, 9.5);

    // announce start of macro
    std::cout << "\n  Beginning cluster-track projection matching macro!" << std::endl;

    // open file w/ frame reader
    podio::ROOTReader reader = podio::ROOTReader();
    reader.openFile( opt.in_file );
    std::cout << "    Opened ROOT-based frame reader." << std::endl;
        
    // get no. of frames and annoucne
    const uint64_t nFrames = reader.getEntries(podio::Category::Event);
    std::cout << "    Starting frame loop: " << reader.getEntries(podio::Category::Event) << " frames to process." << std::endl;

    // iterate through frames (i.e. events in this case)
    uint64_t nClustTotal   = 0;
    uint64_t nClustMatched = 0;
    for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {
        // announce progress
        if(opt.debug) std::cout << "      Processing frame " << iFrame + 1 << "/" << nFrames << "..." << std::endl;

        // grab frame
        auto frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );

        // grab collections
        auto& clusters = frame.get<edm4eic::ClusterCollection>( opt.clusters );
        auto& segments = frame.get<edm4eic::TrackSegmentCollection>( opt.projections );
        auto& mc_particles = frame.get<edm4hep::MCParticleCollection>( "MCParticles" );
        NClustersPerEvent->Fill(clusters.size());
        if (opt.debug)  std::cout << "      Number of MCParticles in frame: " << mc_particles.size() << std::endl;
        
        // loop over MCParticle
        TLorentzVector Particle;
        for (edm4hep::MCParticle mc_part : mc_particles) {
            if (mc_part.getGeneratorStatus() != 1) continue;

            int pdg = mc_part.getPDG();              
            auto momentum = mc_part.getMomentum();   
            
            Particle.SetPxPyPzE(momentum.x,momentum.y,momentum.z,mc_part.getEnergy());
            
            MCMomentum->Fill(Particle.P());
            MCEta->Fill(Particle.Eta());
            MCPt->Fill(Particle.Pt());
            MCPhi->Fill(Particle.Phi());


            if (opt.debug) {
                std::cout << "      [MC Particle] PDG: " << pdg 
                        << " | E: " << Particle.E() << " GeV"
                        << " | Eta: " << Particle.Eta() 
                        << " | Phi: " << Particle.Phi() << std::endl;
            }
        }

        std::map<int32_t, int> trackMatchCount;

        // loop over clusters
        for (size_t iClust = 0; edm4eic::Cluster cluster : clusters) {

            int tracksInsideCut = 0;

            // grab eta/phi of cluster
            const double etaClust  = edm4hep::utils::eta( cluster.getPosition() );
            const double phiClust  = edm4hep::utils::angleAzimuthal( cluster.getPosition() );

            // match based on eta/phi dstiance
            double distMatch = std::numeric_limits<double>::max(); 
            double dEta = std::numeric_limits<double>::max(); 
            double dPhi = std::numeric_limits<double>::max(); 


            // loop over projections to find matching one
            std::optional<edm4eic::TrackPoint> match;
            std::optional<edm4eic::TrackSegment> matchedSegment;
            int SegID=0,BestSegID=-1;
            for (edm4eic::TrackSegment segment : segments) {
                for (edm4eic::TrackPoint projection : segment.getPoints()) {

                // ignore if not pointing to calo or at face of calo
                const bool isInSystem = (projection.system == opt.system);
                const bool isAtFace   = (projection.surface == 1);
                if (!isInSystem || !isAtFace) continue;

                // grab eta/phi of projection
                const double etaProject = edm4hep::utils::eta( projection.position );
                const double phiProject = edm4hep::utils::angleAzimuthal( projection.position );

                

                // get distance to projection
                const double distProject = std::hypot(
                    etaClust - etaProject,
                    DeltaPhi(phiClust, phiProject)
                );

                // if smallest distance found, update variables accordingly
                if (distProject < distMatch) {
                    distMatch = distProject;
                    dEta = etaClust - etaProject;
                    dPhi = DeltaPhi(phiClust,phiProject);
                    match = projection;
                    matchedSegment = segment;
                    BestSegID=SegID;

                }

                if (distProject < opt.MaxREtaPhi) {
                        tracksInsideCut++;
                }

                }  // end point loop
            SegID++;
            }  // end segment loop 

            NTracksMatchedToCluster->Fill(tracksInsideCut);

            if( distMatch<opt.MaxREtaPhi) DeltaEtaPhiHist->Fill(dEta,dPhi);
            REtaPhiHist->Fill(distMatch);
            ++nClustTotal;

            // do analysis if match was found...
            if ( match.has_value() && distMatch<opt.MaxREtaPhi ) {
                edm4eic::TrackPoint matchedProject = match.value();
                if(opt.debug) std::cout << "        Found match! match = " << matchedProject << std::endl;
                auto TrackMomentum = matchedProject.momentum;
                double pTrack = std::sqrt(TrackMomentum.x * TrackMomentum.x + 
                              TrackMomentum.y * TrackMomentum.y + 
                              TrackMomentum.z * TrackMomentum.z);
                
                //High energy limit pTrack=eTrack
                const double EnergyRatio=pTrack/cluster.getEnergy();
                EClusteroverETrack->Fill(EnergyRatio);
                EMCoverECluster->Fill(Particle.E()/cluster.getEnergy());
                EMCoverETrack->Fill(Particle.E()/pTrack);



                const double pT_Track = std::hypot(TrackMomentum.x,TrackMomentum.y);
                EnergyRatiovsTrackPt->Fill(pT_Track,EnergyRatio);
                
                const double etaProject = edm4hep::utils::eta( matchedProject.position );
                EnergyRatiovsTrackEta->Fill(etaProject,EnergyRatio);

                ++nClustMatched;
                
                trackMatchCount[BestSegID]++;

                MCMomentumHited->Fill(Particle.P());
                MCEtaHited->Fill(Particle.Eta());
                MCPtHited->Fill(Particle.Pt());
                MCPhiHited->Fill(Particle.Phi());
                break;
            }

        }  // end cluster loop

        for (const auto& [trackID, count] : trackMatchCount)  if(trackID!=-1) NClustersMatchedToTrack->Fill(count);

    }  // end frame loop
    std::cout << "    Finished frame loop: " << nClustMatched << "/" << nClustTotal << " clusters matched." << std::endl;

    TProfile* EnergyRatioProfilevsTrackPt  = EnergyRatiovsTrackPt->ProfileX();
    EnergyRatioProfilevsTrackPt->SetName("EnergyRatioProfilevsTrackPt");
    EnergyRatioProfilevsTrackPt->SetTitle("Mean E_{Cl}/E_{Tr} vs Track p_{T};p_{T} [GeV];<E_{Cl}/E_{Tr}>");

    TProfile* EnergyRatioProfilevsTrackEta  = EnergyRatiovsTrackEta->ProfileX();
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


    Output->Close();
    
    // announce end and exit
    std::cout << "  End of macro!\n" << std::endl;
    return;

}

// end ------------------------------------------------------------------------