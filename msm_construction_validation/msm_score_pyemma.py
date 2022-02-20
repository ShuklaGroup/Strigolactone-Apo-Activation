import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyemma
import glob
import pickle

def load_data(feat_identifiers):

    file_lists = []
    for ft in feat_identifiers:
        file_lists.append(sorted(glob.glob("../analysis/%s"%ft)))

    feat_all = []
    for i in range(len(file_lists[0])): #Iterate through features
        feat = [] #Initialize feature for j-th trajectory
        for j in range(len(file_lists)):
            feat.append(np.load(file_lists[j][i])) #Load i-th feature for j-th trajectory
        feat_all.append(np.hstack(feat))

    return feat_all

def single_run(feat_all, n_tica_components=10, n_clusters=100):

    tica_object = pyemma.coordinates.tica(data=feat_all, lag=4, dim=n_tica_components)
    tica_trajs = tica_object.get_output()

    cluster_object = pyemma.coordinates.cluster_mini_batch_kmeans(tica_trajs, max_iter=300, k=n_clusters)
    dtrajs = cluster_object.assign(tica_trajs)

    msm_object = pyemma.msm.estimate_markov_model(dtrajs, 428, dt_traj='70 ps', score_method='VAMP1', score_k=10)

    score = msm_object.score_cv(dtrajs, score_method='VAMP1', score_k=5)

    return score

def parameter_search(feat_all, tica_list, clusters_list):

    score_results = []
    for n_tica in tica_list:
        for n_cls in clusters_list:
            score = single_run(feat_all, n_tica_components=n_tica, n_clusters=n_cls)
            score_results.append((n_tica, n_cls, score))

    return score_results

if __name__=="__main__":

    #Read in data, combine features into single matrix
    feat_identifiers = []
    feat_identifiers.extend(["*T1-T3-distance-1*", "*T1-T3-distance-2*", "*T1-T3-distance-3*"]) #T1-T3 distances
    feat_identifiers.extend(["*T2-hel-contact-1*", "*T2-hel-contact-2*", "*T2-hel-contact-3*", "*T2-hel-contact-4*"]) #T2 helicity contacts
    feat_identifiers.extend(["*T1-T2-hinge-contact-1*", "*T1-T2-hinge-contact-2*", "*T1-T2-hinge-contact-3*"]) #T1-T2 hinge
    feat_identifiers.extend(["*T1-T2-distance-1*", "*T1-T2-distance-2*", "*T1-T2-distance-3*"]) #T1-T2 distances
    feat_identifiers.extend(["*ftr-loop-dist*1*", "*ftr-loop-dist*2*", "*ftr-loop-dist*3*", "*ftr-loop-dist*4*", "*ftr-loop-dist*5*", "*ftr-loop-dist*6*", "*ftr-loop-dist*7*"]) #Loop distances

    feat_all = load_data(feat_identifiers)
    #score_results = parameter_search(feat_all, [1,2,3,4,5,6,7,8,9,10], [100, 150, 200, 250, 300, 350, 400, 450])
    score_results = parameter_search(feat_all, range(2,11,2), range(100,601,25))
    #score_results = parameter_search(feat_all, [3,4], [100,200])
    pickle.dump(score_results, open("score_results.pkl",'wb'))

