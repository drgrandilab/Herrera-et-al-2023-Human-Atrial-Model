Matlab code of the Herrera et al. (Am J Physiol Heart Circ Physiol. 2023) human atrial model.

This model was developed by extending the Ni et al. human atrial model (Br J Pharmacol, 2020,
available for download on this website) to integrate an updated formulation of SK channel current
based on recordings in human right-atrial cardiomyocytes from nSR and cAF patients (Heijman et al.
Circ Res. 2023).

_____________________________________________________________________________________________________
Contents:

Herrera_2023_main					loads inital conditions and runs the simulation
Herrera_2023_model_SK_expression			excitation-contraction coupling model
function_beat_analysis					calculates biomarkers from model
calcCurrents						extracts currents from model
yfin_AF_1Hz_1SK						cAF inital conditions (100% GSK, 1 Hz)
yfin_nSR_1Hz_1SK					nSR inital conditions (100% GSK, 1 Hz)

Folder "Population-level analysis":

Figure_2_Sensitivity_Analysis_Biomarkers      		plots results shown in Figure 2
Figure_3_Linear_regression_Alternans     		plots results shown in Figure 3C
Figure_4_Linear_and_Logisitic_regression_DADs      	plots results shown in Figure 4C and D
function_logisitic_regression				performs logisitic regression analysis

.mat files in sub-folder "Population-level analysis/Conditions" contain data extracted from populations
of nSR and cAF models
_____________________________________________________________________________________________________

Reference:

N.T. Herrera, X. Zhang, H. Ni, M.M. Maleckar, J. Heijman, D. Dobrev, E. Grandi, S. Morotti. 
Dual effects of the small-conductance Ca2+-activated K+ current on human atrial 
electrophysiology and Ca2+-driven arrhythmogenesis: an in silico study. 
Am J Physiol Heart Circ Physiol. 2023.

Please cite the above paper when using this model.
