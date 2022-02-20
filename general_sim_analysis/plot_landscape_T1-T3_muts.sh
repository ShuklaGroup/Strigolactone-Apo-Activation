#!/bin/bash

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-G121 --ax_labels T1-T3\ distance\ \(nm\) H247-G121\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_T1T3_dloop-H247-G121-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-R173 --ax_labels T1-T3\ distance\ \(nm\) H247-R173\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename AtD14_T1T3_dloop-H247-R173-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy H247-S168 --ax_labels T1-T3\ distance\ \(nm\) H247-S168\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename AtD14_T1T3_dloop-H247-S168-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy G121-Q214 --ax_labels T1-T3\ distance\ \(nm\) G121-Q124\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_T1T3_dloop-G121-Q124-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy G121-S97 --ax_labels T1-T3\ distance\ \(nm\) G121-S97\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_T1T3_dloop-G121-S97-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy R173-D167 --ax_labels T1-T3\ distance\ \(nm\) R173-D167\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_T1T3_dloop-R173-D167-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files T1-T3-distance.npy R173-P169 --ax_labels T1-T3\ distance\ \(nm\) R173-P169\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_T1T3_dloop-R173-P169-weighted.png --weights msm_use/msm_weights.pkl
