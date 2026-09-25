// interface/PhiStarUtils.h

#ifndef analysis_PhiStarAnalyzer_PhiStarUtils_h
#define analysis_PhiStarAnalyzer_PhiStarUtils_h

#include <algorithm>
#include <cmath>

#include "analysis/PhiStarAnalyzer/interface/AnalysisConfig.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace PhiStarUtils {

    double computePhiStar(const reco::Candidate::LorentzVector& e1, const reco::Candidate::LorentzVector& e2);

    bool passCuts(const reco::Candidate::LorentzVector& e1, const reco::Candidate::LorentzVector& e2, const reco::Candidate::LorentzVector& Z, const CutSet& cuts);

    reco::Candidate::LorentzVector smear(const reco::GenParticle& electron, CLHEP::HepRandomEngine& engine);

    double computeHT(std::vector<const reco::GenJet*>& jets);

    reco::Candidate::LorentzVector sumJets(std::vector<const reco::GenJet*>& jets);

    const reco::Candidate* findOrigin(const reco::GenParticle* p);

    std::vector<const reco::GenJet*> cleanJets(const std::vector<reco::GenJet>& jets, std::vector<const reco::GenParticle*>& electrons, double maxEta, double minPt);

}

#endif
