// -*- C++ -*-
//
// Package:    analysis/PhiStarAnalyzer
// Class:      PhiStarPlotter
//
/**\class PhiStarPlotter PhiStarPlotter.cc analysis/PhiStarAnalyzer/plugins/PhiStarPlotter.cc

 Description: Makes plots for PhiStarAnalyzer

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rand Bovard
//         Created:  Mon, 3 Aug 2026 15:36:41 GMT
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"
#include "TSystem.h"

#include "analysis/PhiStarAnalyzer/interface/PhiStarPlotter.h"

PhiStarPlotter::PhiStarPlotter(const std::string& outputDir)
    : outputDir_(outputDir)
{

    std::vector<double> phistarEdges = {
      0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.025, 0.030, 0.040, 0.050,
      0.060, 0.080, 0.100, 0.130, 0.170, 0.220, 0.280, 0.350, 0.440, 0.560,
      0.700, 0.900, 1.150, 1.500, 2.000, 2.600, 3.400, 4.500, 6.000, 8.000};

    h_phistar_ = new TH1D(
        "h_phistar", "#phi*;#phi*;Events",
        (int)phistarEdges.size() - 1, &phistarEdges[0]);

    h_mZ_ = new TH1D("h_mZ", "Z mass;m_{Z};Events", 150, 50.0, 200.0);
    h_qT_ = new TH1D("h_qT", "Z pT;q_{T};Events", 50, 0.0, 50.0);
    h_yZ_ = new TH1D("h_yZ", "Z rapidity;y_{Z};Events", 100, -5.0, 5.0);

    h_e1_pT_ = new TH1D("h_e1_pT", "e1 pT;p_{T};Events", 150, 0, 150);
    h_e1_eta_ = new TH1D("h_e1_eta", "e1 eta;#eta;Events", 100, -5.0, 5.0);
    h_e1_phi_ = new TH1D("h_e1_phi", "e1 phi;#phi;Events", 64, -3.1416, 3.1416);

    h_e2_pT_ = new TH1D("h_e2_pT", "e2 pT;p_{T};Events", 150, 0, 150);
    h_e2_eta_ = new TH1D("h_e2_eta", "e2 eta;#eta;Events", 100, -5.0, 5.0);
    h_e2_phi_ = new TH1D("h_e2_phi", "e2 phi;#phi;Events", 64, -3.1416, 3.1416);

    h_HT_ = new TH1D("h_HT", "HT;HT;Events", 150, 0.0, 300.0);

    h_nJets_ = new TH1D("h_nJets", "Number of Jets;Number of Jets;Events", 51, -0.5, 50.5);

    h_jet_pts_ = new TH1D("h_jet_pts", "All jet pts", 250, 0, 250);

    h_dR_jets_ = new TH1D("h_dR_jets", "Delta R between electrons and closest jets", 100, 0, 1);

    //h_mpt_ = new TH1D("h_mpt", "Missing transverse momentum", 150, 0, 150);

    //res_pte1_ = new TH2D("res_pte1", "e1 pt difference", 
    //        50, 25, 150,
    //        50, -0.04, 0.04);
    //res_pte2_ = new TH2D("res_pte2", "e2 pt difference", 
    //        50, 15, 150,
    //        50, -0.04, 0.04);
    //res_qt_ = new TH2D("res_qt", "qt difference", 
    //        50, 0, 50,
    //        50, -0.1, 0.1);
    //res_phistar_ = new TH2D("res_phistar", "phistar difference", 
    //        50, 0, 50,
    //        50, -0.01, 0.01);

    h_HT_vs_qT_ = new TH2D("h_HT_vs_qT", "HT vs qT;HT;qT", 
                    20, 0.0, 500.0,
                    20, 0.0, 100.0);

    h_vecHT_vs_qT_ = new TH2D("h_vecHT_vs_qT", "vector HT vs qT;vector HT;qT", 
                    20, 0.0, 100.0,
                    20, 0.0, 100.0);

}

void PhiStarPlotter::fill(
        double phistar, 
        const reco::Candidate::LorentzVector Z, 
        const reco::Candidate::LorentzVector e1, 
        const reco::Candidate::LorentzVector e2, 
        double HT,
        reco::Candidate::LorentzVector sum_jets,
        std::vector<const reco::GenJet*> jets, 
        double res_pte1,
        double res_pte2,
        double res_qt,
        double res_phistar,
        std::vector<double> deltas,
        double weight) 
    {

    h_phistar_->Fill(phistar, weight);

    h_mZ_->Fill(Z.mass(), weight);
    h_qT_->Fill(Z.pt(), weight);
    h_yZ_->Fill(Z.Rapidity(), weight);

    h_e1_pT_->Fill(e1.pt(), weight);
    h_e1_eta_->Fill(e1.eta(), weight);
    h_e1_phi_->Fill(e1.phi(), weight);

    h_e2_pT_->Fill(e2.pt(), weight);
    h_e2_eta_->Fill(e2.eta(), weight);
    h_e2_phi_->Fill(e2.phi(), weight);

    h_HT_->Fill(HT, weight);

    h_nJets_->Fill(jets.size(), weight);

    for (const auto& jet : jets) {
        h_jet_pts_->Fill(jet->pt(), weight);
    }

    for (const auto& dR : deltas) {
        h_dR_jets_->Fill(dR, weight);
    }

    //double mpt = sum_jets.pt();
    //h_mpt_->Fill(mpt, weight);

    //res_pte1_->Fill(e1.pt(), res_pte1, weight);
    //res_pte2_->Fill(e2.pt(), res_pte2, weight);
    //res_qt_->Fill(Z.pt(), res_qt, weight);
    //res_phistar_->Fill(Z.pt(), res_phistar, weight);

    h_HT_vs_qT_->Fill(HT, Z.pt(), weight);

    h_vecHT_vs_qT_->Fill(sum_jets.pt(), Z.pt(), weight);

}

void PhiStarPlotter::write() {
    gSystem->mkdir(outputDir_.c_str(), true);
    TFile fout((outputDir_ + "/plots.root").c_str(), "RECREATE");

    h_phistar_->Write();

    h_mZ_->Write();
    h_qT_->Write();
    h_yZ_->Write();

    h_e1_pT_->Write();
    h_e1_eta_->Write();
    h_e1_phi_->Write();

    h_e2_pT_->Write();
    h_e2_eta_->Write();
    h_e2_phi_->Write();

    h_HT_->Write();

    h_nJets_->Write();

    h_jet_pts_->Write();

    h_dR_jets_->Write();

    //h_mpt_->Write();

    //res_pte1_->Write();
    //res_pte2_->Write();
    //res_qt_->Write();
    //res_phistar_->Write();

    h_HT_vs_qT_->Write();

    h_vecHT_vs_qT_->Write();

    fout.Close();
}
