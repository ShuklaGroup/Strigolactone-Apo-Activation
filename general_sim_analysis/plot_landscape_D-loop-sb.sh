#!/bin/bash

python -m plot_free_energy --ax_files salt*distance-1 helical*unfolded --ax_labels D218-K242\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 2 0 1 --savename AtD14_saltbridge1_weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files salt*distance-2 helical*unfolded --ax_labels K217-E244\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 2 0 1 --savename AtD14_saltbridge2_weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files salt*distance-3 helical*unfolded --ax_labels K217-D167\ distance\ \(nm\) T2\ helical\ content --ax_lims 0 2 0 1 --savename AtD14_saltbridge3_weighted.png --weights msm_use/msm_weights.pkl

