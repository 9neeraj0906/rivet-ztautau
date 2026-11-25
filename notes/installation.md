# ---------- On the host PC (Linux, Docker already installed) ----------

# 1. Create a working folder and unpack the project skeleton into it
mkdir ~/riv
mv ~/Downloads/rivet-ztautau.zip ~/riv/
cd ~/riv
unzip rivet-ztautau.zip

# 2. Create a persistent container named "rivet" from the Rivet tutorial image
#    (Rivet + Pythia8 + MadGraph). The image downloads automatically (~1.6 GB) on first run.
#    -v $HOME/riv:/riv     -> host folder ~/riv appears inside the container as /riv
#    -w /riv               -> container opens in /riv
#    -u $(id -u):$(id -g)  -> run as my user so output files aren't root-owned
docker run -it --name rivet -v $HOME/riv:/riv -w /riv -u $(id -u):$(id -g) hepstore/rivet-tutorial:4.1.2-amd64

# 3. Re-enter the container later / leave it
docker start -ai rivet
exit

# ---------- Inside the container: verify the installation ----------
rivet --version                                     # -> rivet v4.1.2
pythia8-main144 -h                                  # -> main144 usage message
rivet --list-analyses | grep ATLAS_2019_I1768911    # -> ATLAS_2019_I1768911  Z pT and Z phi* at 13 TeV
