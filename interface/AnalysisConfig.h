#ifndef analysis_PhiStarAnalyzer_AnalysisConfig_h
#define analysis_PhiStarAnalyzer_AnalysisConfig_h

#include <string>
#include <vector>

struct CutSet {
    double minPte1;
    double minPte2;
    double maxEta;
    double minEta;
    double maxJetEta;
    double minHT;
    double maxqT;
    bool notInGap;
    unsigned int region;
    bool extra;
};

struct AnalysisConfig {
    std::string label;
    std::string outputDir;
    CutSet cuts;
    std::string source = "Gen";

    std::string pdfWeightId;
    int pdfWeightIndex;
    double sumOfWeights;
    bool useNominalOnly;
};

#endif
