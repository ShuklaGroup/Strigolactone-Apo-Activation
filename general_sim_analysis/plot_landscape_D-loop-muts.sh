#!/bin/bash

python -m plot_free_energy --ax_files A216-backbone-distance-1 loop-distance-3 --ax_labels A216-backbone\ distance\ \(nm\) D218-H247\ distance\ \(nm\) --ax_lims 0 2 0 2 --savename AtD14_dloop-mut-old1_weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files salt*distance-3 loop-distance-3 --ax_labels K217-D167\ distance\ \(nm\) D218-H247\ distance\ \(nm\) --ax_lims 0 2 0 2 --savename AtD14_dloop-mut-old2_weighted.png --weights msm_use/msm_weights.pkl

python -m plot_free_energy --ax_files S220-T215 loop-distance-3 --ax_labels S220-T215\ distance\ \(nm\) D218-H247\ \(nm\) --ax_lims 0 2 0 2 --savename AtD14_dloop-mut-old3_weighted.png --weights msm_use/msm_weights.pkl

