import glob
import numpy as np
import mdtraj as md

#targettopfile="../prmtop_stripped/stripped.active1_noh.prmtop"
#targettopfile="stripped.4ih4_withlig_noh.prmtop"
#targettopfile="stripped.striga-CLIM-inactive.prmtop"
targettopfile="stripped.4ih4_modelled_active.prmtop"
#targettopfile="stripped.ShHTL7_withGR24.prmtop"
#targettopfile="stripped.AtD14_withMeCLA.prmtop"
#targettopfile="stripped.striga-DOH-covH.prmtop"
#targettopfile="stripped.AtD14_active_withGR24.prmtop"
#targettopfile="stripped.5CBK.prmtop"
#targettopfile="stripped.ShHTL7_apo_inactive.prmtop"

pairs = [[144,148],[145,149],[146,150]] #AtD14
#pairs = [[146,150],[147,151],[148,152]] #ShHTL7
#pairs = [[147,151],[148,152],[149,153]] #ShHTL5
#pairs = [[135,156],[138,153],[142,149]] #AtD14
#pairs = [[137,158],[140,155],[144,151]] #ShHTL7
labels = ["1","2","3"]

#for file in glob.glob('./stripped/*inactive*round15*xtc'):
for file in glob.glob('./*xtc'):
    t = md.load(file, top=targettopfile)
    #d = md.compute_contacts(t,[[140,155]], scheme='ca')[0]     	# distance between residue pairs on helix T1 and T2 respectively (ShHTL7)
    for j in range(len(labels)):

        d = md.compute_contacts(t, [pairs[j]], scheme='ca')[0]     	# distance between residue pairs on helix T1 and T2 respectively
        n_frames = t.n_frames

        dis = np.empty([n_frames, 1])

        for i in range(n_frames):
          dis[i,0:1]=d[i][0]

        np.save('./analysis/'+file.split('/')[-1]+'_ftr-helix-T1-T2-hinge-contact-%s.npy'%labels[j], dis)
