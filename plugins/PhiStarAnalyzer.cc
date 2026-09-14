// -*- C++ -*-
//
// Package:    analysis/PhiStarAnalyzer
// Class:      PhiStarAnalyzer
//
/**\class PhiStarAnalyzer PhiStarAnalyzer.cc analysis/PhiStarAnalyzer/plugins/PhiStarAnalyzer.cc

 Description: Main event loop module for calculating phistar and plotting relevant variables from MC Z->ee samples

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rand Bovard
//         Created:  Tue, 28 Jul 2026 10:40:37 GMT
//
//

#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "analysis/PhiStarAnalyzer/interface/AnalysisConfig.h"
#include "analysis/PhiStarAnalyzer/interface/RwlWeightMap.h"
#include "analysis/PhiStarAnalyzer/interface/PhiStarUtils.h"
#include "analysis/PhiStarAnalyzer/interface/PhiStarPlotter.h"

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

class PhiStarAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit PhiStarAnalyzer(const edm::ParameterSet&);

    void endJob() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    //edm::EDGetTokenT<reco::GsfElectronCollection> recoElectronsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

    double getEventWeight(const edm::Event&, const AnalysisConfig& cfg) const;

    std::vector<AnalysisConfig> configs_;
    std::vector<std::unique_ptr<PhiStarPlotter>> plotters_;
};

PhiStarAnalyzer::PhiStarAnalyzer(const edm::ParameterSet& iConfig) : 
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    //recoElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("recoElectrons"))),
    genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("ak4GenJetsNoNu"))),
    lheEventToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
    genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")))
{
    auto variations = iConfig.getParameter<std::vector<edm::ParameterSet>>("variations");

    for (const auto& pset : variations) { 
        AnalysisConfig cfg;

        cfg.label = pset.getParameter<std::string>("label");
        cfg.outputDir = pset.getParameter<std::string>("outputDir");
        cfg.source = pset.getParameter<std::string>("electronSource");

        cfg.cuts.minPte1 = pset.getParameter<double>("minPte1");
        cfg.cuts.minPte2 = pset.getParameter<double>("minPte2");
        cfg.cuts.maxEta = pset.getParameter<double>("maxEta");
        cfg.cuts.minEta = pset.getParameter<double>("minEta");
        cfg.cuts.maxJetEta = pset.getParameter<double>("maxJetEta");
        cfg.cuts.minHT = pset.getParameter<double>("minHT");
        cfg.cuts.maxqT = pset.getParameter<double>("maxqT");
        cfg.cuts.notInGap = pset.getParameter<bool>("notInGap");
        cfg.cuts.region = pset.getParameter<unsigned int>("region");
        cfg.cuts.extra = pset.getParameter<bool>("extra");

        cfg.pdfWeightId = pset.getParameter<std::string>("pdfWeightId");
        cfg.pdfWeightIndex = -1;
        cfg.sumOfWeights = 0.0;
        cfg.useNominalOnly = pset.getParameter<bool>("useNominalOnly");

        plotters_.push_back(std::make_unique<PhiStarPlotter>(cfg.outputDir));

        if (!cfg.useNominalOnly) {
            auto it = kWeightIdToIndex.find(cfg.pdfWeightId);
            if (it == kWeightIdToIndex.end()) {
                throw cms::Exception("Configuration")
                    << "pdfWeightId '" << cfg.pdfWeightId
                    << "' not found in RwlWeightMap.h - check cfi.py against pwg-rwl.dat.";
            }
            cfg.pdfWeightIndex = it->second;
        }

        configs_.push_back(cfg);
    }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
}


double PhiStarAnalyzer::getEventWeight(const edm::Event& iEvent, const AnalysisConfig& cfg) const {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(genInfoToken_, genInfo);
    double nominalWeight = genInfo->weight();

    if (cfg.useNominalOnly || cfg.pdfWeightIndex < 0) {
        return nominalWeight;
    }

    edm::Handle<LHEEventProduct> lheInfo;
    iEvent.getByToken(lheEventToken_, lheInfo);
    const auto& weights = lheInfo->weights();

    if (cfg.pdfWeightIndex >= (int)weights.size()) {
        throw cms::Exception("LogicError") << "PDF weight index out of range.";
    }

    return weights[cfg.pdfWeightIndex].wgt;
}


void PhiStarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    //edm::Handle<reco::GsfElectronCollection> recoElectrons;
    //iEvent.getByToken(recoElectronsToken_, recoElectrons);

    edm::Handle<std::vector<reco::GenJet>> genJets;
    iEvent.getByToken(genJetsToken_, genJets);

    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(iEvent.streamID());

    std::vector<const reco::GenParticle*> genElectrons;
    for (const auto& p : *genParticles) {
        if (std::abs(p.pdgId()) == 11 && p.status() == 1) {
            genElectrons.push_back(&p);
        }
    }

    std::vector<const reco::GenParticle*> zElectrons;
    for (const auto& e: genElectrons) {
        const reco::Candidate* origin = PhiStarUtils::findOrigin(e);
        if (origin && origin->pdgId() == 23) {
            zElectrons.push_back(e);
        }
    }

    if (zElectrons.size() != 2) return;

    std::sort(zElectrons.begin(), zElectrons.end(),
    [](const reco::GenParticle* a, const reco::GenParticle* b) {
        return a->pt() > b->pt();
    });

    const reco::GenParticle* e1 = zElectrons[0];
    const reco::GenParticle* e2 = zElectrons[1];

    for (size_t i = 0; i < configs_.size(); ++i) {
        auto& cfg = configs_[i];

        reco::Candidate::LorentzVector p1, p2;
        if (cfg.source == "Gen") {
            p1 = e1->p4(); p2 = e2->p4();
        } else if (cfg.source == "Smeared") {
            p1 = PhiStarUtils::smear(*e1, engine); p2 = PhiStarUtils::smear(*e2, engine);   
        } else if (cfg.source == "Dressed") {
            p1 = e1->p4(); p2 = e2->p4();
            for (auto& genPh : *genParticles) {
                if (genPh.status() != 1 || genPh.pdgId() != 22) continue;
                if (!genPh.isPromptFinalState()) continue;

                double dR1 = reco::deltaR(*e1, genPh);
                double dR2 = reco::deltaR(*e2, genPh);

                if (dR1 < dR2) {
                    if (dR1 < 0.1) {
                        p1 += genPh.p4();
                    }
                } else {
                    if (dR2 < 0.1) {
                        p2 += genPh.p4();
                    }
                }
            }
        } else {
            throw cms::Exception("Configuration") << "Unknown electronSource: " << cfg.source;
        }

        reco::Candidate::LorentzVector Z = p1 + p2;

        if (cfg.cuts.notInGap) {
            if (std::abs(p1.eta()) > 1.4442 && std::abs(p1.eta()) < 1.566) continue;
            if (std::abs(p2.eta()) > 1.4442 && std::abs(p2.eta()) < 1.566) continue;
        }

        if (PhiStarUtils::passCuts(p1, p2, Z, cfg.cuts)) {

            double eventWeight = getEventWeight(iEvent, cfg);
            cfg.sumOfWeights += eventWeight;

            std::vector<const reco::GenJet*> goodJets = PhiStarUtils::cleanJets(*genJets, zElectrons, cfg.cuts.maxJetEta);
            std::vector<double> deltas;
            double minDR1 = 999.0;
            double minDR2 = 999.0;
            for (auto& jet : *genJets) {
                double dr1 = reco::deltaR(p1, jet);
                double dr2 = reco::deltaR(p2, jet);
                if (dr1 < minDR1) {
                    minDR1 = dr1;
                }
                if (dr2 < minDR2) {
                    minDR2 = dr2;
                }
            }
            deltas.push_back(minDR1);
            deltas.push_back(minDR2);

            double phistar = PhiStarUtils::computePhiStar(p1, p2);

            double HT = PhiStarUtils::computeHT(goodJets);
            reco::Candidate::LorentzVector sum_jets = PhiStarUtils::sumJets(goodJets);

            reco::Candidate::LorentzVector unsm_Z = e1->p4() + e2->p4();
            double unsm_phi = PhiStarUtils::computePhiStar(e1->p4(), e2->p4());

            double res_pte1 = (e1->pt() - p1.pt())/e1->pt();
            double res_pte2 = (e2->pt() - p2.pt())/e2->pt();
            double res_qt = (Z.pt() - unsm_Z.pt())/unsm_Z.pt();
            double res_phistar = (phistar-unsm_phi)/unsm_phi;

            if (HT < cfg.cuts.minHT) continue;

            if (cfg.cuts.region == 1) {
                if (Z.pt() >= sum_jets.pt() - 15.0) continue;
            }

            if (cfg.cuts.region == 2) {
                if (Z.pt() <= sum_jets.pt() + 15.0) continue;
            }

            plotters_[i]->fill(phistar, Z, p1, p2, HT, sum_jets, goodJets, res_pte1, res_pte2, res_qt, res_phistar, deltas, eventWeight);
        }
    }

}

void PhiStarAnalyzer::endJob() {
    for (auto& p : plotters_) p->write();
}


void PhiStarAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    
    edm::ParameterSetDescription desc;
    
    desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("ak4GenJetsNoNu", edm::InputTag("ak4GenJetsNoNu"));
    desc.add<edm::InputTag>("lheInfo", edm::InputTag("externalLHEProducer"));
    desc.add<edm::InputTag>("genEventInfo", edm::InputTag("generator"));

    edm::ParameterSetDescription varDesc;
    varDesc.add<std::string>("label", "");
    varDesc.add<std::string>("outputDir", "");
    varDesc.add<std::string>("electronSource", "Gen");

    varDesc.add<double>("minPte1", 0.0);
    varDesc.add<double>("minPte2", 0.0);
    varDesc.add<double>("maxEta", 2.5);
    varDesc.add<double>("minEta", 0.0);
    varDesc.add<double>("maxJetEta", 5.0);
    varDesc.add<double>("minHT", 0.0);
    varDesc.add<double>("maxqT", 1000.0);
    varDesc.add<bool>("notInGap", false);
    varDesc.add<unsigned int>("region", 0);
    varDesc.add<bool>("extra", false);

    varDesc.add<std::string>("pdfWeightId", "");
    varDesc.add<bool>("useNominalOnly", true);

    desc.addVPSet("variations", varDesc);

    descriptions.add("PhiStarAnalyzer", desc);
}

DEFINE_FWK_MODULE(PhiStarAnalyzer);
