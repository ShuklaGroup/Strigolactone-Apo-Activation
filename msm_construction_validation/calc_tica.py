import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import natsort
import pyemma
import glob
import pickle
import time

def load_data(feat_id, path):

    file_lists = []

    for ft in feat_id:
        file_lists.append(natsort.natsorted(glob.glob("%s/analysis/%s"%(path, ft))))

    print(file_lists[0])

    feat_all = []
    for i in range(len(file_lists[0])): #Iterate through features
        feat = [] #Initialize feature for j-th trajectory
        for j in range(len(file_lists)):
            feat.append(np.load(file_lists[j][i])) #Load i-th feature for j-th trajectory
        feat_all.append(np.hstack(feat))
    
    return feat_all

def transform_feat(feat_all, mean):
    
    feat_transformed = []
    for feat in feat_all:
        feat_len = np.shape(feat)[0]
        feat_transformed.append(feat - np.matmul(np.ones((feat_len,1)), mean.reshape((1,20))))

    return feat_transformed

def compute_tica(feat_all, tica_components, save=False, path=None):

    tica_object = pyemma.coordinates.tica(data=feat_all, lag=4, dim=tica_components)
    tica_trajs = tica_object.get_output()

    if save:
        pickle.dump(tica_object, open("%s/tica_analysis/tica_obj_calc.pkl"%path,'wb'))
        pickle.dump(tica_trajs, open("%s/tica_analysis/tica_trajs_calc.pkl"%path,'wb'))

    return tica_object, tica_trajs

def fit_new_data(tica_object, new_data, save=False, path=None):

    new_tica = tica_object.transform(new_data)

    if save:
        #pickle.dump(new_tica, open("%s/tica_analysis/tica_trajs_fitted.pkl"%path, 'wb'))
        pickle.dump(new_tica, open("./tica_trajs_fitted.pkl", 'wb'))

    return new_tica

def find_correlated_feat(corr_matrix, feat_list, cutoff=0.7):
    
    corr_dim = np.shape(corr_matrix)
    feat_list = np.array(feat_list)
    corr_features = []
    corr = []

    for i in range(corr_dim[1]):
        corr_dim = corr_matrix[:,i]
        corr_features.append(feat_list[np.abs(corr_dim) > cutoff])
        corr.append(corr_dim[np.abs(corr_dim) > cutoff])

    return corr_features, corr

if __name__=="__main__":

    #Feature ids list
    feat_identifiers = []
    feat_identifiers.extend(["*T1-T3-distance-1*", "*T1-T3-distance-2*", "*T1-T3-distance-3*"]) #T1-T3 distances
    feat_identifiers.extend(["*T2-hel-contact-1*", "*T2-hel-contact-2*", "*T2-hel-contact-3*", "*T2-hel-contact-4*"]) #T2 helicity contacts
    feat_identifiers.extend(["*T1-T2-hinge-contact-1*", "*T1-T2-hinge-contact-2*", "*T1-T2-hinge-contact-3*"]) #T1-T2 hinge
    feat_identifiers.extend(["*T1-T2-distance-1*", "*T1-T2-distance-2*", "*T1-T2-distance-3*"]) #T1-T2 distances
    feat_identifiers.extend(["*ftr-loop-dist*1*", "*ftr-loop-dist*2*", "*ftr-loop-dist*3*", "*ftr-loop-dist*4*", "*ftr-loop-dist*5*", "*ftr-loop-dist*6*", "*ftr-loop-dist*7*"]) #Loop distances

    D14_path = "/home/jiming/Storage/sda11/AtD14_apo_unbiased"
    HTL7_path = "/home/jiming/JFC-10TB/bdea6630-3a42-4ebb-88cd-13bce784d34d1/ShHTL7-conf-change"
    HTL7_path = "/home/jiming/Storage/ShHTL7_TMD"

    D14_data = load_data(feat_identifiers, D14_path)
    HTL7_data = load_data(feat_identifiers, HTL7_path)

    D14_tica_obj, D14_tica_trajs = compute_tica(HTL7_data, 6, save=True, path=HTL7_path)
    HTL7_tica_obj, HTL7_tica_trajs = compute_tica(HTL7_data, 8, save=False, path=HTL7_path)
    pickle.dump(HTL7_tica_obj, open("HTL7_tica_obj.pkl", 'wb'))
    pickle.dump(HTL7_tica_trajs, open("HTL7_tica_trajs.pkl", 'wb'))
    pickle.dump(D14_tica_obj, open("HTL_tica_obj.pkl", 'wb'))
    pickle.dump(D14_tica_trajs, open("HTL_tica_trajs.pkl", 'wb'))

    D14_fit = fit_new_data(D14_tica_obj, D14_data, save=True, path=D14_path)
    D14_fit = fit_new_data(HTL7_tica_obj, D14_data, save=True, path=D14_path)
    HTL7_fit = fit_new_data(D14_tica_obj, HTL7_data, save=True, path=HTL7_path)
    
    D14_transformed = transform_feat(D14_data, D14_tica_obj.mean)
    print(D14_transformed)
    print(D14_tica_obj.feature_TIC_correlation)
    print(np.matmul(D14_transformed[0], D14_tica_obj.feature_TIC_correlation))
    print(D14_tica_trajs[0])

    D14_corr_feat, D14_corr = find_correlated_feat(D14_tica_obj.feature_TIC_correlation, feat_identifiers)
    print(D14_tica_obj.feature_TIC_correlation)
    np.savetxt("%s/tica_analysis/tica_corr"%D14_path, D14_tica_obj.feature_TIC_correlation)
    print(D14_corr_feat)
    print(D14_corr)
    HTL7_corr_feat, HTL7_corr = find_correlated_feat(HTL7_tica_obj.feature_TIC_correlation, feat_identifiers)
    print(HTL7_tica_obj.feature_TIC_correlation)
    np.savetxt("%s/tica_analysis/tica_corr"%HTL7_path, HTL7_tica_obj.feature_TIC_correlation)
    print(HTL7_corr_feat)
    print(HTL7_corr)
