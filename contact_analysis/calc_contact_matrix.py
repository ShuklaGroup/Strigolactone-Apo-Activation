import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import glob
import pickle

plt.rc('savefig', dpi=300)
matplotlib.rc('font',family='Helvetica-Normal',size=16)

def compute_contact_distances(traj_list, top):

    dist_all = []
    for traj_file in traj_list:
        traj = md.load(traj_file, top=top) #Subsample because memory :(
        distances, contact_list = md.compute_contacts(traj, contacts='all', scheme='closest-heavy', ignore_nonprotein=False)
        dist_all.append(distances)

    dist_all = np.vstack(dist_all)
    
    return contact_list, dist_all

def compute_contact_probs(contact_distances, cutoff=0.45, weights=None):

    #Identify contacts
    contacts = contact_distances.copy() #Copy array
    print(contacts)
    contacts[contact_distances > cutoff] = 0.0 #Not in contact
    contacts[contact_distances < cutoff] = 1.0 #In contact

    print(contacts)
    #Compute ligand contact probability by residue
    if weights is None:
        contact_probs = np.sum(contacts, axis=0)/np.shape(contacts)[0]
    else:
        contact_probs = np.matmul(contacts.T, weights)
        #contact_probs = probs/np.sum(probs)

    return contact_probs

def compute_eq_contact_probs(contact_distances, dtrajs, msm, cutoff=0.45):

    contacts = contact_distances.copy()
    contacts[contact_distances > cutoff] = 0.0 #Not in contact
    contacts[contact_distances < cutoff] = 1.0 #In contact

    eq_dist = msm.pi
    dtrajs = np.hstack(dtrajs)

    contact_prob_per_state = np.zeros((np.shape(contacts)[1], np.shape(eq_dist)[0]))

    print(np.shape(contact_prob_per_state))

    for i in range(np.shape(contacts)[0]): #Count per-residue in each MSM state
        contact_prob_per_state[:,dtrajs[i]] += contacts[i,:].T

    for i in range(np.shape(eq_dist)[0]): #Normalize by MSM state count
        contact_prob_per_state[:,i] /= len(dtrajs[dtrajs==i])

    eq_contact_probs = np.matmul(contact_prob_per_state, eq_dist)

    return eq_contact_probs

def plot_contact_probs(contact_list, contact_probs, savename):

    #x_coords = np.zeros(np.shape(contact_list)[0])
    #y_coords = np.zeros(np.shape(contact_list)[0])

    #for i in range(np.shape(contact_list)[0]):
    #    x_coords[i] = contact_list[i,0]
    #    y_coords[i] = contact_list[i,1]

    n_res = np.max(contact_list+1)
    matrix = np.zeros((n_res, n_res))

    print(np.shape(contact_list))
    print(np.shape(contact_probs))
    for i in range(np.shape(contact_list)[0]):
        matrix[contact_list[i,0], contact_list[i,1]] = contact_probs[i]

    print(matrix)

    fig, ax = plt.subplots()
    #plt.scatter(x_coords, y_coords, s=5, c=contact_probs)
    plt.pcolormesh(range(n_res), range(n_res), matrix)
    plt.xlabel("Residue i")
    plt.ylabel("Residue j")
    ax.set_aspect(1)
    plt.tight_layout()
    plt.savefig("%s.png"%savename)

if __name__=="__main__":

    #contact_list, contact_distances = compute_contact_distances(sorted(glob.glob("../kmc_run/trajs_for_contact_matrix/*xtc")), "../stripped.ShHTL7_apo_active.prmtop")
    #contact_list, contact_distances = compute_contact_distances(sorted(glob.glob("../kmc_run/trajs_for_contact_matrix/*xtc")), "../stripped.4ih4_modelled_active.prmtop")
    #contact_list, contact_distances = compute_contact_distances(["../inactive.xtc"], "../stripped.ShHTL7_apo_active.prmtop")
    #contact_list, contact_distances = compute_contact_distances(sorted(glob.glob("../stripped/*")), "../stripped.4ih4_modelled_active.prmtop")
    #contact_list, contact_distances = compute_contact_distances(["../active.xtc"], "../stripped.4ih4_modelled_active.prmtop")
    #for i in range(400):
    #    contact_list, contact_distances = compute_contact_distances(["../msm_use/state_%d.xtc"%i], "../stripped.4ih4_modelled_active.prmtop")
        #contact_list, contact_distances = compute_contact_distances(["../active.xtc"], "../stripped.4ih4_modelled_active.prmtop")

    #    if i == 0:
    #        np.save("contact_list.npy", contact_list)

    #    np.save("contact_distances_%d.npy"%i, contact_distances)

    #    contact_probs = compute_contact_probs(contact_distances)
    #    np.save("contact_probs_%d.npy"%i, contact_probs)

    contact_list = np.load("contact_list.npy")
    contact_probs = np.load("contact_probs.npy")
    plot_contact_probs(contact_list, contact_probs, "contact_prob-cmesh.png")
