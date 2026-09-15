#ifndef analysis_PhiStarAnalyzer_PhiStarPlotter_h
#define analysis_PhiStarAnalyzer_PhiStarPlotter_h

#include <string>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

class PhiStarPlotter {
public:
    PhiStarPlotter(const std::string& outputDir);
    void fill(
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
            double mpt,
            double weight
            );
    void write();

private:
    std::string outputDir_;
    TH1D* h_phistar_;
    TH1D* h_mZ_;
    TH1D* h_yZ_;
    TH1D* h_qT_;
    TH1D* h_e1_pT_;    
    TH1D* h_e1_eta_;    
    TH1D* h_e1_phi_;    
    TH1D* h_e2_pT_;    
    TH1D* h_e2_eta_;    
    TH1D* h_e2_phi_;    
    TH1D* h_HT_;
    TH1D* h_nJets_;
    TH1D* h_jet_pts_;
    TH1D* h_mpt_;
    TH1D* h_dR_jets_;

    //TH2D* res_pte1_;
    //TH2D* res_pte2_;
    //TH2D* res_qt_;
    //TH2D* res_phistar_;

    TH2D* h_HT_vs_qT_;
    TH2D* h_vecHT_vs_qT_;
};

#endif
