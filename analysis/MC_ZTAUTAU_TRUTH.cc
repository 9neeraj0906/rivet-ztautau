// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {

  /// Truth-level Z->tautau, mu + tau_h channel.
  /// Kinematic cuts follow my thesis (CMS Open Data) mu-tau_h selection,
  /// translated to particle level. Detector-only cuts are not applied.
  class MC_ZTAUTAU_TRUTH : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_ZTAUTAU_TRUTH);

    void init() {

      // Projections

      // Muons: all stable muons, including those from tau decays
      FinalState fsm(Cuts::abseta < 2.1 && Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON);
      declare(fsm, "Muons");

      // Hadronic taus: no cut here, the pT/eta cuts go on the visible tau later
      TauFinder hadronic_taus(TauDecay::HADRONIC);
      declare(hadronic_taus, "HadronicTaus");

      // Electrons for the extra-electron veto (prompt = not from hadron decays,
      // the truth-level stand-in for isolation)
      PromptFinalState fse(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON);
      declare(fse, "Electrons");

      // Jets for the b-jet veto: anti-kT R=0.4, no muons or neutrinos
      FinalState fsj(Cuts::abseta < 4.9);
      FastJets jets(fsj, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jets, "Jets");


      /// Book Histograms

      book(_h["mvis"],       "mvis",       30,  0., 150.);
      book(_h["h_tau_full"], "h_tau_full", 30,  0., 150.);
      book(_h["mu_pt"],      "mu_pt",      20, 20., 100.);
      book(_h["tauh_pt"],    "tauh_pt",    20, 32., 132.);
      book(_h["dR"],         "dR",         25,  0.,   5.);

      /// Fiducial cross-section counter
      book(_c["sigma_fid"], "sigma_fid");
    }

    void analyze(const Event& event) {
      Particles leptons;
      const FinalState& p1 = apply<FinalState>(event, "Muons");
      leptons = p1.particlesByPt();

      Particles Taus;
      const TauFinder& p2 = apply<TauFinder>(event, "HadronicTaus");
      Taus = p2.particlesByPt();

      if (leptons.size() == 1 && Taus.size() == 1 && leptons[0].charge() != Taus[0].charge()) {

        // Visible hadronic tau: sum of its stable decay products, without neutrinos
        FourMomentum vis;
        for (const Particle& d : Taus[0].stableDescendants()) {
          if (isNeutrino(d)) continue;
          vis += d.mom();
        }

        // Thesis tau_h cuts and mu-tau separation
        if (vis.pT() < 32*GeV || vis.abseta() > 2.3) vetoEvent;
        if (deltaR(leptons[0].mom(), vis) < 0.5) vetoEvent;

        // Extra-electron veto: any prompt electron separated from the tau_h by dR > 0.5
        const Particles electrons = apply<PromptFinalState>(event, "Electrons").particles();
        for (const Particle& e : electrons) {
          if (deltaR(e.mom(), vis) > 0.5) vetoEvent;
        }

        // b-jet veto: jets with pT > 20 GeV, |eta| < 2.4 (thesis cuts),
        // b-tagged by a ghost-associated b-hadron with pT > 5 GeV
        const FastJets& p3 = apply<FastJets>(event, "Jets");
        const Jets jets = p3.jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);
        for (const Jet& j : jets) {
          if (j.bTagged(Cuts::pT > 5*GeV)) vetoEvent;
        }

        // ---- Event passed all cuts: fill ----

        // Muon + full hadronic tau (neutrinos included)
        FourMomentum momentum = leptons[0].mom() + Taus[0].mom();
        double h_tau_full = momentum.mass();
        _h["h_tau_full"]->fill(h_tau_full/GeV);

        // Muon + visible hadronic tau
        FourMomentum visMomentum = leptons[0].mom() + vis;
        double mvis = visMomentum.mass();
        _h["mvis"]->fill(mvis/GeV);

        _h["mu_pt"]->fill(leptons[0].pT()/GeV);
        _h["tauh_pt"]->fill(vis.pT()/GeV);
        _h["dR"]->fill(deltaR(leptons[0].mom(), vis));

        _c["sigma_fid"]->fill();
      }

    }

    void finalize() {
      const double sf = crossSection()/picobarn/sumW();
      for (auto hist : _h) {scale(hist.second, sf);}
      for (auto cnt  : _c) {scale(cnt.second,  sf);}
    }

    /// Histograms and counters
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;

  };

  RIVET_DECLARE_PLUGIN(MC_ZTAUTAU_TRUTH);

}
