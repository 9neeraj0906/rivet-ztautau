# Notes: from zero to a truth-level Z → ττ Rivet analysis

These notes record everything I did, in order, including what went wrong and how I fixed it. They're written so that another student starting with Rivet can follow them as a tutorial.

**Starting point:** comfortable with Python, C++ (Geant4), ROOT, uproot/awkward; M.Sc. thesis on Z → ττ with CMS Open Data (NanoAOD). Never used Rivet, YODA, HepMC, Pythia or MadGraph before.

---

## Contents
1. [The big picture](#1-the-big-picture)
2. [Environment setup (Docker)](#2-environment-setup-docker)
3. [First events with Pythia8](#3-first-events-with-pythia8)
4. [Learning Rivet: the LHC tutorial](#4-learning-rivet-the-lhc-tutorial)
5. [Getting a Z → ττ sample](#5-getting-a-z--ττ-sample)
6. [Writing my own routine](#6-writing-my-own-routine)
7. [Results and physics checks](#7-results-and-physics-checks)
8. [Problems I hit and how I fixed them](#8-problems-i-hit-and-how-i-fixed-them)
9. [What I'd do next](#9-what-id-do-next)
10. [Resources](#10-resources)

---

## 1. The big picture

```
MadGraph / Pythia8  ──►  HepMC events  ──►  Rivet + my routine  ──►  YODA histograms  ──►  plots
 "generate collisions"    "event record"     "select, fill"           "results"           rivet-mkhtml
```

- **HepMC**: a generator-independent text record of every particle in every event, with mother/daughter links. Truth (particle) level, no detector.
- **Rivet**: runs analysis routines on HepMC events. Published LHC measurements are preserved as Rivet routines, so any generator can be compared with data. Uses **projections** (reusable, cached observable calculators like `FinalState`, `FastJets`, `TauFinder`).
- **YODA**: Rivet's histogram format (plain text).
- A Rivet routine is a C++ class with three methods: `init()` (declare projections, book histograms), `analyze()` (once per event: apply projections, cut, fill), `finalize()` (normalise).

Why Rivet matters: experiments unfold measurements to particle level inside a *fiducial region*. Rivet preserves the exact selection, so theorists can compare predictions, and tools like **Contur** can check whether a new-physics model would already have shown up in existing measurements.

---

## 2. Environment setup (Docker)

Image: `hepstore/rivet-tutorial:4.1.2-amd64` (Rivet 4.1.2, Pythia8, MadGraph5_aMC@NLO 3.6.7).

```bash
mkdir ~/riv
docker run -it --name rivet -v $HOME/riv:/riv -w /riv -u $(id -u):$(id -g) hepstore/rivet-tutorial:4.1.2-amd64
# later:
docker start -ai rivet
```
- `-v $HOME/riv:/riv`: my PC folder appears inside the container, so all files live on my disk.
- `-u $(id -u):$(id -g)`: output files are owned by me, not root.
- For installing packages: `docker exec -u root -it rivet bash` from the host.

Check: `rivet --version` → `rivet v4.1.2`; `pythia8-main144 -h` prints usage.

---

## 3. First events with Pythia8

A Pythia8 command file for Drell–Yan Z → ee/μμ (`generation/pythia8/zll_13TeV.cmnd`):

```
Beams:eCM = 13000
WeakSingleBoson:ffbar2gmZ = on
PhaseSpace:mHatMin = 50.
PhaseSpace:mHatMax = 150.
23:onMode = off
23:onIfAny = 11 13
Main:writeHepMC = on
Next:numberCount = 100
```
```bash
pythia8-main144 -c zll_13TeV.cmnd -n 1000 -o outputs/zll_test
```
- The mass window cuts away the low-mass γ* region (which rises like 1/m²) and is wider than any analysis window, so final-state radiation doesn't bias the edges.
- **Lesson:** no spaces around the colon. `PhaseSpace: mHatMin` is silently not recognised.
- Sanity check on the first event: x₁ = 0.106, x₂ = 4.7×10⁻⁴, so √(x₁x₂s) = 91.7 GeV, matching the event's scale. The event record agrees with basic kinematics.
- 1000 full LHC events ≈ 190 MB of HepMC text (~930 particles per event).

---

## 4. Learning Rivet: the LHC tutorial

I followed the Rivet LHC tutorial (C. Gütschow, Z+jets truth analysis). The version bundled in the Docker image (`/work/TUTORIAL.cc`) is already written for Rivet 4. My solutions are in `tutorial/`.

| Exercise | Concept | What I saw |
|---|---|---|
| 1 | `FinalState` + `Cuts`, analysis options (`LMODE`) | Z peak at 91 GeV; electron curve depleted at the peak vs muons |
| 2 | Dressed leptons (`LeptonFinder`, ΔR = 0.1) | Dressing recovers FSR photons; e and μ curves agree much better |
| 3 | Jets (`FastJets`, anti-kT R = 0.4), multiplicity, H_T | Electron channel has extra "jets": electrons get clustered, muons are excluded |
| 4 | Overlap removal (`discardIfAnyDeltaRLess`) | After removing jets near leptons, e and μ channels agree |
| 5 | b-tagging by ghost association (`hasBTag`) | b-jet multiplicity |

Exercise 6 (scale-variation uncertainty bands) I skipped for now; I'll revisit it after learning more MC generator basics.

**Rivet 3 → 4 differences I hit** (the PDF version of the tutorial is Rivet 3):
| Rivet 3 | Rivet 4 |
|---|---|
| `DEFAULT_RIVET_ANALYSIS_CTOR` | `RIVET_DEFAULT_ANALYSIS_CTOR` |
| `DECLARE_RIVET_PLUGIN` | `RIVET_DECLARE_PLUGIN` |
| `sumOfWeights()` | `sumW()` |
| `DressedLeptons` | `LeptonFinder` (photon/lepton argument order swapped) |
| `ZFinder` | `DileptonFinder` |
| `FastJets::ANTIKT`, `JetAlg::Muons::NONE` | `JetAlg::ANTIKT`, `JetMuons::NONE` |
| `rivet-mkhtml --errs --no-weights` | these flags no longer exist |

Tip: `rivet-mkanalysis NAME` generates a Rivet 4 template; its comments are a good syntax reference.

---

## 5. Getting a Z → ττ sample

- The tutorial's Z+jets sample has **no** Z → ττ (checked by counting PDG ±15 particles: 0).
- **ATLAS Open Data** released HepMC event-generation samples in 2025. With `atlasopenmagic` (release `2025r-evgen-13tev`, 5173 datasets) I found MadGraph+Pythia8 Z → ττ samples (`..._Ztautau_Np0`…`Np4`, e.g. DSID 361510). The https download links returned **HTTP 404**; I reported this to ATLAS Open Data support.
- So I generated my own sample with MadGraph + Pythia8. See section 8 for the problems along the way.

**Final sample:** 10,000 events of pp → Z/γ* → τ⁺τ⁻ at 13 TeV, leading order (MadGraph 3.6.7, NNPDF2.3 LO), m(ττ) > 50 GeV, showered and hadronised by Pythia8, which also decays the taus. MadGraph cross-section 1493 ± 1.5 pb; Rivet's `MC_XS` reads 1491.7 pb (agreement within 0.1%). Cards are in `generation/madgraph/`.

---

## 6. Writing my own routine

Goal: translate my thesis μτ_h selection to particle level.

**Which thesis cuts carry over, and which don't:**
| Thesis (reco level) | Truth level |
|---|---|
| μ pT > 20 GeV, \|η\| < 2.1 | same, on the truth muon |
| τ_h pT > 32 GeV, \|η\| < 2.3 | same, on the **visible** τ_h |
| opposite sign, ΔR(μ,τ_h) > 0.5 | same |
| extra-electron veto | prompt electrons pT > 10 GeV, \|η\| < 2.5, ΔR > 0.5 from τ_h |
| b-jet veto (DeepCSV) | anti-kT 0.4 jets pT > 20, \|η\| < 2.4, ghost-associated b-hadron pT > 5 GeV |
| MET filters, PV, trigger + matching, IDs, isolation, dz/dxy, pileup | no truth equivalent: dropped (prompt requirement stands in for isolation) |

**Key physics point: the visible tau.** `TauFinder` returns the full tau, *including* its neutrino. A detector only sees the visible decay products, so:
```cpp
FourMomentum vis;
for (const Particle& d : tau.stableDescendants()) {
  if (isNeutrino(d)) continue;
  vis += d.mom();
}
```
The τ_h cuts and the visible mass use `vis`, not the full tau.

**Development stages** (each one checked with a plot before moving on):
- **A:** m(μ + full τ_h), no τ_h cuts. Peaks ~70–75 GeV, below m_Z, because the leptonic tau's neutrinos are still missing.
- **B:** added m_vis(μ + visible τ_h). Peak moves to ~60 GeV with a long low-mass tail: neutrinos on both legs.
- **C:** thesis τ_h cuts + opposite sign + ΔR. Low-mass tail removed.
- **Final:** + electron veto, b-jet veto, more histograms, σ_fid counter.

**Fiducial cross-section:** a Rivet counter (`CounterPtr`) filled once per selected event and scaled by σ_gen/ΣW gives σ_fid directly in pb. Acceptance A = σ_fid / σ_gen.

---

## 7. Results and physics checks

- **m_vis** peaks at about 60–70 GeV, well below 91 GeV, and is broad, because neutrinos from both tau decays escape with a different energy in every event.
- **Selection sculpting:** the τ_h pT > 32 GeV and μ pT > 20 GeV cuts remove low-mass events, so the surviving distribution sits higher than without cuts.
- **ΔR(μ, τ_h)** peaks near π: taus from a Z decay are roughly back-to-back.
- **Muon pT** falls steeply from 20 GeV: the muon shares its tau's energy with two neutrinos.
- **Acceptance cross-check:** only 2 × BR(τ→μνν) × BR(τ→hadrons) ≈ 2 × 0.174 × 0.648 ≈ 0.23 of Z → ττ events are μτ_h, so the acceptance must be below 0.23. Mine: **(0.72 ± 0.09)%**, from 72 selected events, well below the 0.23 μτ_h ceiling.

---

## 8. Problems I hit and how I fixed them

### Docker image environment
- **PyROOT segfaults on `import ROOT`.** This also broke `import yoda` and `rivet-mkhtml`. Reinstalling `g++-11` (ROOT's expected compiler) removed the first error but not the crash; running as root didn't help either. Workaround: a dummy module so the optional ROOT import fails cleanly:
  ```bash
  mkdir -p /riv/noroot
  echo 'raise ImportError("ROOT disabled")' > /riv/noroot/ROOT.py
  export PYTHONPATH=/riv/noroot:$PYTHONPATH
  ```
- **`rivet-mkhtml` LaTeX errors** (`type1cm.sty`, `type1ec.sty` not found): installed `texlive-latex-extra` and `cm-super` as root, then ran `texhash`. Check with `kpsewhich type1cm.sty`.

### Rivet / C++ mistakes worth knowing
- **Read only the first compiler error**; the long "candidate: …" lists are noise.
- `"Found double-booking of /…"`: two histograms booked with the same output name.
- **Segmentation fault:** indexing past the end of a list (`Taus[1]` when there's one tau). Check `.size()` first.
- **Always re-run `rivet-build` after editing**; otherwise Rivet runs the old `.so`.
- Library names must match `Rivet*.so`. The class name must match the analysis name used with `-a`.
- An editor tab pointing at a file you've `mv`ed will silently save a new copy in the old place.
- `rivet -n 100` reads only 100 events. Don't plot a test run by mistake.
- Projections can't select taus with `FinalState`: taus are unstable. Use `TauFinder`.

### MadGraph + Pythia8
- **Reweighting module fails to compile** (`FileNotFoundError: 'meson'`): on Python ≥ 3.12, numpy's `f2py` needs `meson`, which isn't in the image. Left `reweight` OFF (not needed).
- **Shower step fails with exit code 127** (`main164: symbol lookup error: undefined symbol …Pythia8…Settings…`): the container has two Pythia8 installations: system Pythia 8.315 in `/usr/local`, and MadGraph's copy in `HEPTools/pythia8`. MadGraph compiles `main164` against its own copy, but `LD_LIBRARY_PATH` made it load the system `libpythia8.so`. Diagnosed with `ldd`. Fix: an alias that sets the library path only for MadGraph:
  ```bash
  alias mg5='LD_LIBRARY_PATH=/work/MG5_aMC_v3_6_7/HEPTools/pythia8/lib:$LD_LIBRARY_PATH /work/MG5_aMC/bin/mg5_aMC'
  ```
  A session started before the alias still fails. I re-showered the existing run instead of regenerating it (`launch ztautau -i`, then `pythia8 run_03`).
- **Default generator-level lepton cuts apply to taus too** (`ptl`, `etal`, `drll` in `run_card.dat`). Undecayed taus don't correspond to any analysis cut, so this can bias the sample (cross-section 843.9 pb with defaults). Set `ptl = 0`, `etal = -1`, `drll = 0`, kept `mmll = 50`. **Rule:** every event that could pass the analysis cuts must also pass the generator cuts.
- **Interactive menu timers** (60/90 s) can skip your `set` commands. Edit `Cards/run_card.dat` in an editor before `launch`.
- Output: HepMC2 (`IO_GenEvent`, 2.06.09), read directly by Rivet 4. 146 weights per event (nominal + 145 scale/PDF variations); `--skip-weights` runs on the nominal only. ~670 MB gzipped for 10k events.

### ATLAS Open Data
- `atlasopenmagic dataset urls … --protocol https` returns URLs prefixed with `simplecache::` (strip it). They gave HTTP 404 on both `opendata.cern.ch` and `eospublic.cern.ch` (the latter also needs CERN's Grid CA certificate). Reported to atlas-outreach-opendata-support@cern.ch.

---

## 9. What I'd do next
- Larger sample (50–100k events) for smoother distributions.
- The both-full-taus mass as a check (should peak at m_Z).
- Other channels (eτ_h, τ_hτ_h) via an analysis option, and a check of the tau decay-mode fractions (≈35% leptonic, ≈65% hadronic).
- Compare with an NLO or merged sample (e.g. ATLAS Open Data, once downloadable).
- Tutorial Exercise 6: scale-variation uncertainty bands from the event weights.
- Extend to HH → bb̄ττ, and try Contur.

---

## 10. Resources

**Rivet**
- Rivet homepage: https://heprivet.gitlab.io/
- Rivet tutorials (incl. "What is a Projection?", "Writing a basic analysis", "Migrating from v3 to v4"): https://heprivet.gitlab.io/tutorials/
- LHC basic tutorial (Z+jets, dressing, jets, b-tags, weights): https://heprivet.gitlab.io/tutorials/lhc-basic-tutorial
- Doxygen code reference: https://heprivet.gitlab.io/doxy/latest/
- Analysis list with source code: https://heprivet.gitlab.io/analyses/
- Running Rivet via Docker: https://heprivet.gitlab.io/tutorials/docker/
- Rivet 4 release paper: https://arxiv.org/abs/2404.15984
- Rivet 3 paper: https://arxiv.org/abs/1912.05451

**Docker images**
- `hepstore/rivet-tutorial` (Rivet + Pythia8 + MadGraph): https://hub.docker.com/r/hepstore/rivet-tutorial

**Generators**
- Pythia 8 online manual: https://pythia.org/latest-manual/Welcome.html
- MadGraph5_aMC@NLO: https://launchpad.net/mg5amcnlo
- A. Buckley et al., "General-purpose event generators for LHC physics": https://arxiv.org/abs/1101.2599

**Open data**
- ATLAS Open Data event-generation (HepMC) samples: https://opendata.atlas.cern/docs/data/for_research/evgen_data
- CERN Open Data record for the 13 TeV HEPMC release: https://opendata.cern.ch/record/160000
- "Open LHC Monte Carlo Event Generation" (2026): https://arxiv.org/abs/2605.12229

**Physics reference values**
- Tau branching fractions and particle codes: Particle Data Group, https://pdg.lbl.gov
