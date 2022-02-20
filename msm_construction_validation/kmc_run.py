import numpy as np
import mdtraj as md
import pyemma
import pickle
import argparse
import glob

def run_all(msm, initial_state, nsteps, topfile, traj_list, save_name="kmc_traj"):

    tmatrix = msm.P
    state_traj = kmc(msm.P, initial_state, nsteps)
    print(state_traj)
    np.save("%s.npy"%save_name, state_traj)
    structure_traj = get_structure_traj(state_traj, msm, topfile, sorted(glob.glob("./%s/*"%traj_list)))
    structure_traj.save_xtc("%s.xtc"%save_name)

def kmc(tmatrix, initial_state, nsteps):

    nstates = np.shape(tmatrix)[0]
    moves = np.random.rand(nsteps)
    state_traj = np.zeros(nsteps, dtype=int)
    state_traj[0] = initial_state
    #print(np.shape(tmatrix))
    for i in range(nsteps):
        current_state = state_traj[i]
        #print(current_state)
        j = 0
        while np.sum(tmatrix[current_state,:j+1]) < moves[i]:
            j += 1
        state_traj[i] = j

    return state_traj

def get_structure_traj(state_traj, msm_obj, topfile, traj_list, split=False):

    #Get samples from each MSM state
    all_state_samples = msm_obj.sample_by_state(10000)
    #print(all_state_samples)
    #print(all_state_samples[0])
    #Iterate through state_traj, load random frame from state
    if split:
        n_trajs = int(len(state_traj)/1000)

        for i in range(n_trajs):
            frames = []
            for j in range(1000):
                if np.mod(j,1) == 0:
                    state = state_traj[i*1000+j]
                    print(state)
                    smp = all_state_samples[state][np.random.randint(10000)]
                    #print(len(traj_list))
                    #print(smp[0])
                    #print(smp[1])
                    frames.append(md.load_frame(traj_list[smp[0]], smp[1], top=topfile))

            structure_traj = md.join(frames)
            structure_traj.save_xtc("kmc_traj_%d.xtc"%i)

    else:
        frames = []
        for i in range(len(state_traj)):
            if np.mod(i,1) == 0:
                state = state_traj[i] 
                print(state)
                smp = all_state_samples[state][np.random.randint(10000)] #Draw sample
                frames.append(md.load_frame(traj_list[smp[0]], smp[1], top=topfile))
     
        structure_traj = md.join(frames)
        structure_traj.save_xtc("kmc_traj.xtc")
        
    return structure_traj

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--initial_state", type=int, help="Initial MSM state")
    parser.add_argument("--n_steps", type=int, help="Number of steps to run")
    parser.add_argument("--msm_path", type=str, help="Path to MSM output")
    parser.add_argument("--traj_path", type=str, help="Path to trajectories")
    parser.add_argument("--topfile", type=str, help="Path to topology file")
    args = parser.parse_args()

    return args

if __name__=="__main__":
    
    args = get_args()

    msm = pickle.load(open("%s/msm_object.pkl"%args.msm_path, 'rb'))
    if args.initial_state is None:
        initial_state = np.random.randint(msm.n_states)
        run_all(msm, initial_state, args.n_steps, args.topfile, args.traj_path, save_name="kmc_traj")
    else:
        run_all(msm, args.initial_state, args.n_steps, args.topfile, args.traj_path, save_name="kmc_traj")
    kmc_traj = np.load("kmc_traj.npy")
    get_structure_traj(kmc_traj, msm, args.topfile, sorted(glob.glob("%s/*.xtc"%args.traj_path)), split=True)
