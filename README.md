# Bytnerowicz-Akana-Griffin-Menge-N-fix-Temp-Code
R scripts to recreate analyses and figures for Bytnerowicz, Akana, Griffin, and Menge "The temperature sensitivity of woody dinitrogen fixation across species and growing temperatures"


R scripts are split into two folders: "Figures" and "Analysis"

Data used in analyses and for figures are available in the "Data" folder. R scripts are annotated to indicate what folder required data are stored in.



#For "Figures":

Script for figures found in the main text is Figures_Main_Text.R

Scripts for supplementary figures have the supplementary figure number in the filename (e.g., SFigures_1,3,4_R_Script.R is for supplementary figures 1, 3, and 4)



#For "Analysis":

Scripts to calculate N-fixation rates:

Km_Calculations_Alnus.R, Km_Calculations_Gliricidia.R, Km_Calculations_Morella.R, and Km_Calculations_Robinia.R calculate the Michaelis-Menten half-saturation constant (Km) and 95% CI for measuring nitrogenase activity with ARACAS with (Supplementary Table 4) and without a possible change of Km with measurement temperature (Supplementary Table 5)

CF_Calculations.R calculates the Conversion Factor (CF) values in Supplementary Table 13 and examines the possible effect of temperature on CF values


Scripts to model the temperature response of N-fixation:

Nfix_Model_Comparison.R compares models for N-fixation (Supplementary Tables 6 and 8) and then generates parameter estimates and 95% CIs for N-fixation (Supplementary Tables 1 and 3)

Nfix_LinKm_Estimates.R generates parameter estimates and 95% CIs for N-fixation with the possibility of linear acclimation of Km (Supplementary Tables 11 and 12)


Scripts to calculate and model the temperature response of carbon exchange:

Resp_Calculations.R tests for the effect of growing temperature on leaf respiration in the light (RL), compares acclimation models for the effect of temperature on leaf respiration (Supplementary Table 10), and calculates the relative respiration rate at 25 deg. C for A-Ci curve calculations

A-Ci_Calculations.R fits the A-Ci curves in order to calculate A275, Vcmax, and Jmax

Photosynthesis_Model_Comparison.R compares models for photosynthesis (both A275 and Asat; Supplementary Tables 7 and 9) and then generates parameter estimates and 95% CIs for photosynthesis (Supplementary Table 2)


Scripts to calculate the effect of prolonged exposure of N-fixation to high temperatures:

Piecewise_Reg_Constant_Temp_Fits.R fits piecewise functions (equation S4) to N-fixation data at constant temperatures

Temporal_Decline_EqS7_Parameter_Calculations.R estimates the parameter values in Equation S7

Temporal_Decline_Nfix_Minus_Temporal_Effect.R simulates N-fixation data had they not been affected by the time spent at each temperature and refits the modified beta function (Equation 5) to these simulated data

Nfix_Curves_Prolonged_Exposure.R simulates N-fixation temperature response curves for different heating rates
