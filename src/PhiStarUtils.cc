// -*- C++ -*-
//
// Package:    analysis/PhiStarAnalyzer
// Class:      PhiStarUtils
//
/**\class PhiStarUtils PhiStarUtils.cc analysis/PhiStarAnalyzer/plugins/PhiStarUtils.cc

 Description: functions for PhiStarAnalyzer

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rand Bovard
//         Created:  Mon, 3 Aug 2026 15:58:17 GMT
//
//

#include "analysis/PhiStarAnalyzer/interface/PhiStarUtils.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace PhiStarUtils {

    bool passCuts(
            const reco::Candidate::LorentzVector& e1, 
            const reco::Candidate::LorentzVector& e2, 
            const reco::Candidate::LorentzVector& Z, 
            const CutSet& cuts)
    {
        if (
                e1.pt() > cuts.minPte1 && 
                e2.pt() > cuts.minPte2 && 
                std::abs(e1.eta()) < cuts.maxEta && 
                std::abs(e1.eta()) > cuts.minEta && 
                std::abs(e2.eta()) < cuts.maxEta && 
                std::abs(e2.eta()) > cuts.minEta && 
                Z.pt() < cuts.maxqT) 
        {
            return true;
        }
        return false;
    }


    double computePhiStar(
            const reco::Candidate::LorentzVector& e1, 
            const reco::Candidate::LorentzVector& e2) 
    {
        double dPhi = reco::deltaPhi(e1.phi(), e2.phi());
        double phiAcop = M_PI - std::abs(dPhi);

        double dEta = e1.eta() - e2.eta();
        double sinThetaStar = 1.0 / std::cosh(dEta / 2.0);

        return std::tan(phiAcop / 2.0) * sinThetaStar;
    }


    reco::Candidate::LorentzVector smear(const reco::GenParticle& electron, CLHEP::HepRandomEngine& engine) {

        double a = 0.0;
        double b = 0.0;
        double c = 0.0;

        if (std::abs(electron.eta()) < 1.4442) {
            a = 0.028;
            b = 0.120;
            c = 0.003;
        } else {
            a = 0.057;
            b = 0.770;
            c = 0.0055;
        }

        double E = electron.energy();
        double err = std::sqrt( (a*a)/E + (b*b)/(E*E) + c*c );

        double smearFactor = CLHEP::RandGaussQ::shoot(&engine, 1.0, err);
        reco::Candidate::LorentzVector p4(
                electron.px()*smearFactor,
                electron.py()*smearFactor,
                electron.pz()*smearFactor,
                E*smearFactor
                );

        return p4;
    }


    double computeHT(std::vector<const reco::GenJet*>& jets) {
        double ht = 0.0;
        for (const auto& jet : jets) {
            ht += jet->pt();
        }
        return ht;
    }


    reco::Candidate::LorentzVector sumJets(std::vector<const reco::GenJet*>& jets) {
        reco::Candidate::LorentzVector summedJet(0, 0, 0, 0);
        for (auto& jet : jets) {
            summedJet += jet->p4();
        }
        return summedJet;
    }


    const reco::Candidate* findOrigin(const reco::GenParticle* p) {
        const reco::Candidate* mom = p->mother();
        while (mom && std::abs(mom->pdgId()) == 11) {
            mom = mom->mother();
        }
        return mom;
    }


    std::vector<const reco::GenJet*> cleanJets(
            const std::vector<reco::GenJet>& jets, 
            std::vector<const reco::GenParticle*>& electrons,
            double maxEta) 
    {
        std::vector<const reco::GenJet*> goodJets;
        for (const auto& jet : jets) {
            bool isolated = std::all_of(electrons.begin(), electrons.end(),
                [&jet](const reco::GenParticle* e) {
                    return reco::deltaR(jet, *e) > 0.1;
                });
            if (isolated && std::abs(jet.eta()) < maxEta) {
                goodJets.push_back(&jet);
            }

        }
        return goodJets;
    }
}
