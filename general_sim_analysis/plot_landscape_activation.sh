#!/bin/bash

python -m plot_free_energy --ax_files helix-T1-T3-distance.npy helical-content-unfolded --ax_labels T1-T3\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 4 0 1 --savename AtD14_apo_weighted-1.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files helix-T1-T3-distance.npy T1-T2-hinge-contact-1 --ax_labels T1-T3\ distance\ \(nm\) T1-T2\ hinge\ contact\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_apo_weighted-hinge.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files helix-T1-T3-distance.npy T1-T2-hinge-contact-2 --ax_labels T1-T3\ distance\ \(nm\) N151-K155-distance --ax_lims 0 4 0 2 --savename AtD14_apo_weighted-hinge-2.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files helix-T1-T3-distance.npy T1-T2-hinge-contact-3 --ax_labels T1-T3\ distance\ \(nm\) N152-K156-distance --ax_lims 0 4 0 2 --savename AtD14_apo_weighted-hinge-3.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files helix-T1-T3-distance.npy loop-distance-3 --ax_labels T1-T3\ distance\ \(nm\) D-loop\ distance\ \(nm\) --ax_lims 0 4 0 3 --savename AtD14_apo_weighted-dloop.png --weights msm_use/msm_weights.pkl
