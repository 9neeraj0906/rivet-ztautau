BEGIN PLOT /MC_ZTAUTAU_TRUTH/mvis
Title=Visible mass, $\mu\tau_h$ channel (truth level)
XLabel=$m_\mathrm{vis}(\mu,\tau_h)$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_\mathrm{vis}$ [pb/GeV]
LogY=0
END PLOT

BEGIN PLOT /MC_ZTAUTAU_TRUTH/h_tau_full
Title=Muon + full hadronic tau (neutrinos included)
XLabel=$m(\mu,\tau_h^\mathrm{full})$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m$ [pb/GeV]
LogY=0
END PLOT

BEGIN PLOT /MC_ZTAUTAU_TRUTH/mu_pt
Title=Muon transverse momentum
XLabel=$p_T^\mu$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_T$ [pb/GeV]
LogY=0
END PLOT

BEGIN PLOT /MC_ZTAUTAU_TRUTH/tauh_pt
Title=Visible hadronic tau transverse momentum
XLabel=$p_T^{\tau_h,\mathrm{vis}}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_T$ [pb/GeV]
LogY=0
END PLOT

BEGIN PLOT /MC_ZTAUTAU_TRUTH/dR
Title=Angular separation between muon and hadronic tau
XLabel=$\Delta R(\mu,\tau_h)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta R$ [pb]
LogY=0
END PLOT
