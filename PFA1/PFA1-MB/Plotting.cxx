#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

void SetLargeAxisStyle(TH1* h) {
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.15);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
}

void SetLargeAxisStyle(TH2* h) {
    SetLargeAxisStyle(static_cast<TH1*>(h));
    h->GetZaxis()->SetTitleSize(0.05);
    h->GetZaxis()->SetLabelSize(0.045);
    h->GetZaxis()->SetTitleOffset(1.1);
}

void SetCommonStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(0.045);
}

TH1* GetHist1D(TFile* file, const char* name) {
    TH1* h = dynamic_cast<TH1*>(file->Get(name));
    if (!h) {
        std::cerr << "ERROR: Failed to load histogram '" << name << "' from file '" << file->GetName() << "'\n";
    }
    return h;
}

TH2* GetHist2D(TFile* file, const char* name) {
    TH2* h = dynamic_cast<TH2*>(file->Get(name));
    if (!h) {
        std::cerr << "ERROR: Failed to load histogram '" << name << "' from file '" << file->GetName() << "'\n";
    }
    return h;
}

void DrawSpeciesPlots(const char* label, TFile* file, Int_t color) {
    TH1* MCMomentum = GetHist1D(file, "MCMomentum");
    TH1* MCEta = GetHist1D(file, "MCEta");
    TH1* MCPt = GetHist1D(file, "MCPt");
    TH1* MCPhi = GetHist1D(file, "MCPhi");
    TH1* MCMomentumHited = GetHist1D(file, "MCMomentumHited");
    TH1* MCEtaHited = GetHist1D(file, "MCEtaHited");
    TH1* MCPtHited = GetHist1D(file, "MCPtHited");
    TH1* MCPhiHited = GetHist1D(file, "MCPhiHited");
    TH2* DeltaEtaPhiHist = GetHist2D(file, "DeltaEtaPhiHist");
    TH1* REtaPhiHist = GetHist1D(file, "REtaPhiHist");
    TH1* EClusteroverETrack = GetHist1D(file, "EClusteroverETrack");
    TH1* EMCoverECluster = GetHist1D(file, "EMCoverECluster");
    TH1* EMCoverETrack = GetHist1D(file, "EMCoverETrack");
    TH1* EClustMinusFsubEtrk = GetHist1D(file, "EClustMinusFsubEtrk");
    TH2* EnergyRatiovsTrackPt = GetHist2D(file, "EnergyRatiovsTrackPt");
    TH2* EnergyRatiovsTrackEta = GetHist2D(file, "EnergyRatiovsTrackEta");
    TProfile* EnergyRatioProfilevsTrackPt = dynamic_cast<TProfile*>(file->Get("EnergyRatioProfilevsTrackPt"));
    TProfile* EnergyRatioProfilevsTrackEta = dynamic_cast<TProfile*>(file->Get("EnergyRatioProfilevsTrackEta"));
    TH1* NTracksMatchedToCluster = GetHist1D(file, "NTracksMatchedToCluster");
    TH1* NClustersMatchedToTrack = GetHist1D(file, "NClustersMatchedToTrack");

    if (!MCMomentum || !MCEta || !MCPt || !MCPhi || !MCMomentumHited || !MCEtaHited || !MCPtHited || !MCPhiHited ||
        !DeltaEtaPhiHist || !REtaPhiHist || !EClusteroverETrack || !EMCoverECluster || !EMCoverETrack ||
        !EnergyRatiovsTrackPt || !EnergyRatiovsTrackEta || !EnergyRatioProfilevsTrackPt || !EnergyRatioProfilevsTrackEta ||
        !NTracksMatchedToCluster || !NClustersMatchedToTrack) {
        std::cerr << "WARNING: One or more histograms missing for " << label << ". Skipping plots for this species.\n";
        return;
    }

    MCMomentumHited->SetLineColor(color);
    MCMomentumHited->SetMarkerColor(color);
    MCEtaHited->SetLineColor(color);
    MCEtaHited->SetMarkerColor(color);
    MCPtHited->SetLineColor(color);
    MCPtHited->SetMarkerColor(color);
    MCPhiHited->SetLineColor(color);
    MCPhiHited->SetMarkerColor(color);
    REtaPhiHist->SetLineColor(color);
    REtaPhiHist->SetMarkerColor(color);
    EClusteroverETrack->SetLineColor(color);
    EClusteroverETrack->SetMarkerColor(color);
    EMCoverECluster->SetLineColor(color);
    EMCoverECluster->SetMarkerColor(color);
    EMCoverETrack->SetLineColor(color);
    EMCoverETrack->SetMarkerColor(color);
    EClustMinusFsubEtrk->SetLineColor(color);
    NTracksMatchedToCluster->SetLineColor(color);
    NTracksMatchedToCluster->SetMarkerColor(color);
    NClustersMatchedToTrack->SetLineColor(color);
    NClustersMatchedToTrack->SetMarkerColor(color);
    EnergyRatioProfilevsTrackPt->SetLineColor(color);
    EnergyRatioProfilevsTrackPt->SetMarkerColor(color);
    EnergyRatioProfilevsTrackPt->SetMarkerStyle(20);
    EnergyRatioProfilevsTrackEta->SetLineColor(color);
    EnergyRatioProfilevsTrackEta->SetMarkerColor(color);
    EnergyRatioProfilevsTrackEta->SetMarkerStyle(20);

    TCanvas* c_truth_hited = new TCanvas(Form("c_truth_hited_%s", label), Form("MC truth vs MC hitted - %s", label), 1400, 1100);
    c_truth_hited->Divide(2,2);

    c_truth_hited->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    SetLargeAxisStyle(MCMomentum);
    MCMomentum->SetTitle("MC truth momentum vs MC hitted; p [GeV]; Events");
    MCMomentum->SetLineColor(kBlue);
    MCMomentum->SetMarkerColor(kBlue);
    MCMomentum->Draw("hist");
    MCMomentumHited->SetLineColor(kRed);
    MCMomentumHited->SetMarkerColor(kRed);
    MCMomentumHited->Draw("hist same");
    
    TLegend leg(0.80, 0.91, 0.98, 0.98);
        leg.AddEntry(MCMomentum, "MC truth", "l");
        leg.AddEntry(MCMomentumHited, "MC hitted", "l");
    leg.Draw("same");
    

    c_truth_hited->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    SetLargeAxisStyle(MCEta);
    MCEta->SetTitle("MC truth #eta vs MC hitted #eta;#eta; Events");
    MCEta->SetLineColor(kBlue);
    MCEta->SetMarkerColor(kBlue);
    MCEta->Draw("hist");
    MCEtaHited->SetLineColor(kRed);
    MCEtaHited->SetMarkerColor(kRed);
    MCEtaHited->Draw("hist same");

    leg.Draw("same");
    

    c_truth_hited->cd(3);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    SetLargeAxisStyle(MCPt);
    MCPt->SetTitle("MC truth p_{T} vs MC hitted p_{T};p_{T} [GeV]; Events");
    MCPt->SetLineColor(kBlue);
    MCPt->SetMarkerColor(kBlue);
    MCPt->Draw("hist");
    MCPtHited->SetLineColor(kRed);
    MCPtHited->SetMarkerColor(kRed);
    MCPtHited->Draw("hist same");
    leg.Draw("same");

    c_truth_hited->cd(4);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    SetLargeAxisStyle(MCPhi);
    MCPhi->SetTitle("MC truth #phi vs MC hitted #phi;#phi; Events");
    MCPhi->SetLineColor(kBlue);
    MCPhi->SetMarkerColor(kBlue);
    MCPhi->SetMinimum(0);
    MCPhi->Draw("hist");
    MCPhiHited->SetLineColor(kRed);
    MCPhiHited->SetMarkerColor(kRed);
    MCPhiHited->Draw("hist same");
    leg.Draw("same");
    c_truth_hited->Update();
    c_truth_hited->SaveAs(Form("Plots/%s/c_mc_truth_hited_%s.pdf", label, label));

    TCanvas* c_delta = new TCanvas(Form("c_delta_%s", label), Form("DeltaEtaPhi and R_{EtaPhi} - %s", label), 1500, 700);
    c_delta->Divide(2,1);
    c_delta->cd(1);
    gPad->SetRightMargin(0.13);
    SetLargeAxisStyle(DeltaEtaPhiHist);
    DeltaEtaPhiHist->SetTitle("Cluster-Track matching #Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi");
    DeltaEtaPhiHist->Draw("COLZ");
    c_delta->cd(2);
    SetLargeAxisStyle(REtaPhiHist);
    REtaPhiHist->SetTitle("Matching distance R_{#eta#phi};R_{#eta#phi};Events");
    REtaPhiHist->Draw("hist");
    c_delta->Update();
    c_delta->SaveAs(Form("Plots/%s/c_delta_%s.pdf", label, label));

    TCanvas* c_eOver = new TCanvas(Form("c_eOver_%s", label), Form("Energy ratios - %s", label), 1000, 800);
    c_eOver->Divide(1,3);
    c_eOver->cd(1);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(EClusteroverETrack);
    EClusteroverETrack->SetTitle("Cluster energy / track momentum;E_{Cl}/p_{Tr};Events");
    EClusteroverETrack->Draw("hist");
    c_eOver->cd(2);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(EMCoverECluster);
    EMCoverECluster->SetTitle("MC energy / cluster energy;E_{MC}/E_{Cl};Events");
    EMCoverECluster->Draw("hist");
    c_eOver->cd(3);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(EMCoverETrack);
    EMCoverETrack->SetTitle("MC energy / track momentum;E_{MC}/p_{Tr};Events");
    EMCoverETrack->Draw("hist");
    c_eOver->Update();
    c_eOver->SaveAs(Form("Plots/%s/c_eOver_%s.pdf", label, label));

    TCanvas* c_2d_profiles = new TCanvas(Form("c_2d_profiles_%s", label), Form("Energy ratios with profiles - %s", label), 1500, 700);
    c_2d_profiles->Divide(2,1);
    c_2d_profiles->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(EnergyRatiovsTrackPt);
    EnergyRatiovsTrackPt->SetTitle("E_{Cl}/E_{Tr} vs track p_{T};p_{T} [GeV];E_{Cl}/E_{Tr}");
    EnergyRatiovsTrackPt->Draw("COLZ");
    EnergyRatioProfilevsTrackPt->SetLineColor(kRed);
    EnergyRatioProfilevsTrackPt->SetMarkerColor(kRed);
    EnergyRatioProfilevsTrackPt->SetLineWidth(2);
    EnergyRatioProfilevsTrackPt->Draw("same");
    {
        TLegend leg(0.62, 0.72, 0.90, 0.88);
        leg.SetTextSize(0.04);
        leg.AddEntry(EnergyRatiovsTrackPt, "2D ratio map", "f");
        leg.AddEntry(EnergyRatioProfilevsTrackPt, "Mean profile", "l");
        leg.Draw();
    }
    c_2d_profiles->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    SetLargeAxisStyle(EnergyRatiovsTrackEta);
    EnergyRatiovsTrackEta->SetTitle("E_{Cl}/E_{Tr} vs track #eta;#eta;E_{Cl}/E_{Tr}");
    EnergyRatiovsTrackEta->Draw("COLZ");
    EnergyRatioProfilevsTrackEta->SetLineColor(kRed);
    EnergyRatioProfilevsTrackEta->SetMarkerColor(kRed);
    EnergyRatioProfilevsTrackEta->SetLineWidth(2);
    EnergyRatioProfilevsTrackEta->Draw("same");
    {
        TLegend leg(0.62, 0.72, 0.90, 0.88);
        leg.SetTextSize(0.04);
        leg.AddEntry(EnergyRatiovsTrackEta, "2D ratio map", "f");
        leg.AddEntry(EnergyRatioProfilevsTrackEta, "Mean profile", "l");
        leg.Draw();
    }
    c_2d_profiles->Update();
    c_2d_profiles->SaveAs(Form("Plots/%s/c_2d_profiles_%s.pdf", label, label));

    TCanvas* c_matched = new TCanvas(Form("c_matched_%s", label), Form("Matched multiplicity - %s", label), 1500, 700);
    c_matched->Divide(2,1);
    c_matched->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(NTracksMatchedToCluster);
    NTracksMatchedToCluster->SetTitle("Matched tracks per cluster;N_{tracks};Events");
    NTracksMatchedToCluster->Draw("hist");
    c_matched->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    SetLargeAxisStyle(NClustersMatchedToTrack);
    NClustersMatchedToTrack->SetTitle("Matched clusters per track;N_{clusters};Events");
    NClustersMatchedToTrack->Draw("hist");
    c_matched->Update();
    c_matched->SaveAs(Form("Plots/%s/c_matched_%s.pdf", label, label));

    TCanvas* c_EClustMinusFsubEtrk = new TCanvas(Form("c_EClustMinusFsubEtrk_%s", label), Form("Energy ratios - %s", label), 1000, 800);
    gPad->SetBottomMargin(0.15);
    SetLargeAxisStyle(EClustMinusFsubEtrk);
    EClustMinusFsubEtrk->Rebin(4);
    EClustMinusFsubEtrk->SetTitle("E_{clust} - f_{sub}E_{trk};E_{clust} - f_{sub}E_{trk} [GeV];Events");
    EClustMinusFsubEtrk->Draw("hist");
    c_EClustMinusFsubEtrk->Update();
    c_EClustMinusFsubEtrk->SaveAs(Form("Plots/%s/c_EClustMinusFsubEtrk_%s.pdf", label, label));
}

void Plotting()
{
    SetCommonStyle();

    TFile* ElectronBarrelFile = TFile::Open("Plots/ElectronBarrel.root", "READ");
    TFile* ElectronNEndCapFile = TFile::Open("Plots/ElectronNEndCap.root", "READ");
    TFile* ElectronPEndCapFile = TFile::Open("Plots/ElectronPEndCap.root", "READ");
    TFile* PionBarrelFile = TFile::Open("Plots/PionBarrel.root", "READ");
    TFile* PionNEndCapFile = TFile::Open("Plots/PionNEndCap.root", "READ");
    TFile* PionPEndCapFile = TFile::Open("Plots/PionPEndCap.root", "READ");
    TFile* PionLFHCALFile = TFile::Open("Plots/PionLFHCAL.root", "READ");

    DrawSpeciesPlots("ElectronBarrel", ElectronBarrelFile, kBlue);
    DrawSpeciesPlots("ElectronNEndCap", ElectronNEndCapFile, kBlue);
    DrawSpeciesPlots("ElectronPEndCap", ElectronPEndCapFile, kBlue);
    DrawSpeciesPlots("PionBarrel", PionBarrelFile, kRed);
    DrawSpeciesPlots("PionNEndCap", PionNEndCapFile, kRed);
    DrawSpeciesPlots("PionPEndCap", PionPEndCapFile, kRed);
    DrawSpeciesPlots("PionLFHCAL", PionLFHCALFile, kRed);

}