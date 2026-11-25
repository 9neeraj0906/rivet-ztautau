// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// Tutorial analysis with dilepton and jet observables
  class TUTORIAL : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TUTORIAL);


    /// Book histograms and initialise projections before the run
    void init() {

      // Electron or muon mode? (default to electron)
      _eemode = (getOption("LMODE") != "MU");

      // Projections

      FinalState fse(Cuts:: abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON);
      declare(fse, "Electrons");
      FinalState fsm(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::MUON);
      declare(fsm, "Muons");

      LeptonFinder dressed_eles(Cuts:: abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON, 0.1, LeptonOrigin::NODECAY);
      declare(dressed_eles, "Dressed Electrons");
      LeptonFinder dressed_muons(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::MUON, 0.1, LeptonOrigin::NODECAY);
      declare(dressed_muons, "Dressed Muons");

      const FinalState fsj(Cuts::abseta < 4.9);
      FastJets jets(fsj, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jets, "Jets");


      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll_bare"],    "mll_bare", mll_bins);
      book(_h["mll_dressed"], "mll_dressed", mll_bins);
      book(_h["njet_all"],    "njet_all",   6, -0.5,  5.5);
      book(_h["njet_iso"],    "njet_iso",   6, -0.5,  5.5);
      book(_h["nbjet_all"],   "nbjet_all",  3, -0.5,  2.5);
      book(_h["nbjet_iso"],   "nbjet_iso",  3, -0.5,  2.5);
      book(_h["jet1pt"],      "jet1pt",      6,  20., 110.);
      book(_h["HT"],          "HT",         6,  20., 150.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Get leptons (bare & dressed)
      Particles leptons;
      if (_eemode == true){
        const FinalState& p1 = apply<FinalState>(event, "Electrons");
        leptons = p1.particlesByPt();
      } else{
        const FinalState& p1 = apply<FinalState>(event, "Muons");
        leptons = p1.particlesByPt();
      }

      if (leptons.size() == 2){
        if (leptons[0].charge() != leptons[1].charge()){
          FourMomentum momentum = leptons[0].mom() + leptons[1].mom();
          double invm = momentum.mass();
          _h["mll_bare"]->fill(invm/GeV);
        }
      }
      Particles dressedLeptons;
      if (_eemode == true){
        const LeptonFinder& p1 = apply<LeptonFinder>(event, "Dressed Electrons");
        dressedLeptons = p1.particlesByPt();
      } else{
        const LeptonFinder& p1 = apply<LeptonFinder>(event, "Dressed Muons");
        dressedLeptons = p1.particlesByPt();
      }

      if (dressedLeptons.size() == 2){
        if (dressedLeptons[0].charge() != dressedLeptons[1].charge()){
          FourMomentum momentum = dressedLeptons[0].mom() + dressedLeptons[1].mom();
          double invm_dressed = momentum.mass();
          Jets fastjets;
          if (invm_dressed > 66*GeV && invm_dressed < 116*GeV){
            const FastJets& p1 = apply<FastJets>(event,"Jets");
            fastjets = p1.jetsByPt(Cuts::absrap < 4.5 && Cuts::pT>10*GeV);
            _h["njet_all"]->fill(fastjets.size());
            int bjets_all;
            bjets_all = count(fastjets, hasBTag());
            _h["nbjet_all"]->fill(bjets_all);
            Jets isojets = discardIfAnyDeltaRLess(fastjets, dressedLeptons, 0.4);
            _h["njet_iso"]->fill(isojets.size());
            int bjets_iso;
            bjets_iso = count(isojets, hasBTag());
            _h["nbjet_iso"]->fill(bjets_iso);
            double HT = 0;
            for (const Jet& j:isojets){
              HT = HT + j.pT();
            }
            _h["HT"]->fill(HT/GeV);
          }
          _h["mll_dressed"]->fill(invm_dressed/GeV);
        }
      }
    }


    /// Post-run histogram manipulations
    void finalize() {

      // Normalise histograms
      const double sf = crossSection()/picobarn / sumW();
      for (auto hist : _h) { scale(hist.second, sf); }
    }


    /// Histograms
    map<string, Histo1DPtr> _h;

    /// Mode flag
    bool _eemode;

  };


  RIVET_DECLARE_PLUGIN(TUTORIAL);

}
