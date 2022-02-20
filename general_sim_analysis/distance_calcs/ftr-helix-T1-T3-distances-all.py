import glob
import numpy as np
import mdtraj as md

#targettopfile="stripped.4ih4_modelled_active.prmtop"
targettopfile="stripped.ShHTL7_apo_active.prmtop"

#pairs = [[136,166],[139,171],[143,175]] #AtD14
pairs = [[138,167],[141,172],[145,176]] #ShHTL7
labels = ["1","2","3"]

#for file in glob.glob('./stripped/*inactive*round15*.xtc'):
for file in glob.glob('./*.xtc'):
    t = md.load(file, top=targettopfile)
    #d = md.compute_contacts(t,[[140,155]], scheme='ca')[0]     	# distance between residue pairs on helix T1 and T2 respectively (ShHTL7)
    for j in range(len(labels)):

        d = md.compute_contacts(t, [pairs[j]], scheme='ca')[0]     	# distance between residue pairs on helix T1 and T2 respectively
        n_frames = t.n_frames

        dis = np.empty([n_frames, 1])

        for i in range(n_frames):
          dis[i,0:1]=d[i][0]

        np.save('./analysis/'+file.split('/')[-1]+'_ftr-helix-T1-T3-distance-%s.npy'%labels[j], dis)
