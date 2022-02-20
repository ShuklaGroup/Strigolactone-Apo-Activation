#!/bin/bash

#python -m plot_free_energy --ax_files S220 helical*unfolded --ax_labels S220-T215\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 2 0 1 --savename AtD14_S220_weighted.png --weights msm_use/msm_weights.pkl

#python -m plot_free_energy --ax_files  S220 helical*d-loop --ax_labels S220-T215\ distance\ \(nm\) D-loop\ helical\ content --ax_lims 0 2 0 1 --savename AtD14_S220_dloop-hel-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files  helix-T1-T3-distance.npy helical*d-loop --ax_labels S220-T215\ distance\ \(nm\) D-loop\ helical\ content --ax_lims 0 4 0 1 --savename AtD14_T1T3_dloop-hel-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files  S220-H247 helical*d-loop --ax_labels S220-T215\ distance\ \(nm\) D-loop\ helical\ content --ax_lims 0 4 0 1 --savename AtD14_S220_dloop-hel-weighted.png --weights msm_use/msm_weights.pkl

