from ROOT import *
import numpy as n
import sys
import os.path
from math import ceil

file_id = sys.argv[1]
# the hadd of all the output ntuples
fname = file_id
nick = os.path.basename(fname).split("_")[-1].rsplit(".")[0]
pt = [20, 26, 30, 34]
numberOfHLTTriggers = 6

saveOnlyOS = True # True; save only OS, False: save both and store weight for bkg sub

#######################################################
fIn = TFile.Open(fname)
#tIn = TChain('Ntuplizer/TagAndProbe')
tIn = fIn.Get('Ntuplizer/TagAndProbe')
#tTriggerNames = TChain("Ntuplizer/triggerNamesProbe")
tTriggerNames = fIn.Get("Ntuplizer/triggerNamesProbe")
#files = glob.glob(iname)
puw_file = TFile.Open("/storage/9/mburkart/Ntuplizer_SingleTau_v2/pu_weightfile.root", "read")
puw_histo = puw_file.Get(nick)

outname = fname.replace ('.root', '_forFit.root')
fOut = TFile (outname, 'recreate')
tOut = tIn.CloneTree(0)
tOutNames = tTriggerNames.CloneTree(-1) # copy all

briso   = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
brnoiso = [n.zeros(1, dtype=int) for x in range (0, len(pt))]
bkgSubW = n.zeros(1, dtype=float)
puweight = n.zeros(1, dtype=float)

hltPathTriggered_OS   = [n.zeros(1, dtype=int) for x in range (0, numberOfHLTTriggers+1)]

for i in range (0, len(pt)):
    name = ("hasL1_" + str(pt[i]))
    tOut.Branch(name, brnoiso[i], name+"/I")
    name += "_iso"
    tOut.Branch(name, briso[i], name+"/I")

for i in range (0, numberOfHLTTriggers):
    tTriggerNames.GetEntry(i)
    name = ("hasHLTPath_" + str(i))
    tOut.Branch(name, hltPathTriggered_OS[i], name+"/I")

tOut.Branch("isoHLT", hltPathTriggered_OS[6], name+"/I")

tOut.Branch("bkgSubW", bkgSubW, "bkgSubW/D")
tOut.Branch("puweight", puweight, "puweight/D")

nentries = tIn.GetEntries()
for ev in range (0, nentries):
    tIn.GetEntry(ev)
    if (ev%10000 == 0) : print ev, "/", nentries

    if abs(tIn.tauProbeEta) > 2.1:
        continue

    if saveOnlyOS and not tIn.isOS:
        continue
    
    if tIn.eleTagPt < 40.:
        continue

    if tIn.tauTrkPt < 50.:
        continue

    #if tIn.tauProbeDM != 10:
    #    continue

    if tIn.tauProbePt < 50.:
        continue

    if tIn.isProbeHLTmatched == 0:
        continue

    for i in range (0, len(pt)):
        briso[i][0] = 0
        brnoiso[i][0] = 0

    for i in range (0, numberOfHLTTriggers):
        hltPathTriggered_OS[i][0] = 0

    L1iso = True if tIn.l1tIso == 1 else False
    L1pt = tIn.l1tPt
    for i in range(0, len(pt)):
        # print L1pt, pt[i]
        #
        if L1pt > pt[i]:
            brnoiso[i][0] = 1
            # print "SUCCESS!! ", brnoiso[i]
            if L1iso:
                briso[i][0] = 1

    triggerBits = tIn.eleProbeTriggerBits
    HLTpt = tIn.hltPt
    for bitIndex in range(0, numberOfHLTTriggers):
        if ((triggerBits >> bitIndex) & 1) == 1:
            hltPathTriggered_OS[bitIndex][0] = 1

    bkgSubW[0] = 1. if tIn.isOS else -1.

    if (L1pt > 26) and (L1iso) and (HLTpt > 32) and (((triggerBits >> 2) & 1) == 1):
        hltPathTriggered_OS[6][0] = 1
    else:
        hltPathTriggered_OS[6][0] = 0

    puweight[0] = puw_hist.GetBinContent(int(ceil(tIn.nTruePU))) if (tIn.nTruePU >= 0. and tIn.nTruePU <= 200.) else 0.

    tOut.Fill()

tOutNames.Write()
tOut.Write()
fOut.Close()
fIn.Close()
