#!/bin/bash

python -m plot_free_energy --ax_files loop-distance-3 H247-G121 --ax_labels D218-H247\ distance\ \(nm\) H247-G121\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_dloop-H247-G121-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 H247-R173 --ax_labels D218-H247\ distance\ \(nm\) H247-R173\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename AtD14_dloop-H247-R173-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 H247-S168 --ax_labels D218-H247\ distance\ \(nm\) H247-S168\ distance\ \(nm\) --ax_lims 0 4 0 4 --savename AtD14_dloop-H247-S168-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 G121-Q214 --ax_labels D218-H247\ distance\ \(nm\) G121-Q124\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_dloop-G121-Q124-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 G121-S97 --ax_labels D218-H247\ distance\ \(nm\) G121-S97\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_dloop-G121-S97-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 R173-D167 --ax_labels D218-H247\ distance\ \(nm\) R173-D167\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_dloop-R173-D167-weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files loop-distance-3 R173-P169 --ax_labels D218-H247\ distance\ \(nm\) R173-P169\ distance\ \(nm\) --ax_lims 0 4 0 2 --savename AtD14_dloop-R173-P169-weighted.png --weights msm_use/msm_weights.pkl
