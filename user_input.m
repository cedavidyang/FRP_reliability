% user input file
% constant are defined by user in this file.
clear; clc;

% give the datafile path and name
XLS_DATAFILE_PATH_AND_NAME = '.\shear_database.xlsx';
% XLS_DATAFILE_PATH_AND_NAME = '.\flexure_database.xlsx';                          

% sub-database of interest
% shear+side: shear strengthening with side bonded FRP

SUB_TEST_DATABASE_NAME ='shear+side';
% SUB_TEST_DATABASE_NAME ='flexure+IC';

% design code used
DESIGN_CODE = 'hk';

%% basic information for the establishment of design cases (shear)
% H_DESIGN_MM: height of beam in mm
% AS_DESIGN_MM: distance from steel centroid to concrete surface in mm
% H2B_DESIGN: ratio of beam height to beam width
% FC_DESIGN_MPA: nominal concrete compressive strength in MPa
% BAR_TYPE_DESIGN: type of steel bars, 1 for deformed bar
%                                      2 for round bar
% SD_DESIGN_MM: bar diameter in mm
% SS_MIN_DESIGN_MM: mimimum interval of steel bars in mm
% SS_MAX_DESIGN_MM: maximum interval of steel bars in mm
% S2D_DESIGN: shear span ratio
% FS_DESIGN_MPA: norminal yielding strength of steel in MPa
% FRP_TYPE_DESIGN: type of FRP strengthening, 1 for CFRP
%                                             2 for GFRP
% BETA_DESIGN_DEG: inclination of FRP in degree
% T_FRP_DESIGN_MM: thickness of FRP in mm
% STRENGTH_TYPE: strengthening type, 1 for continuous strengthening
%                                    2 for discrete strengthening
% E_CARBON_DESIGN_MPA: elastic modulus of CFRP in MPa
% STRAIN_CARBON_DESIGN: ultimate strain of CFRP
% E_GLASS_DESIGN_MPA: elastic modulus of GFRP in MPa
% STRAIN_GLASS_DESIGN: ultimate strain of GFRP
% FRP_CONFIG_DESIGN: FRP configuration, 1 for Sheet
%                                       2 for Strip
% FRP_FORM_DESIGN: form of FRP strengthing 1 for Side
%                                          2 for U
%                                          3 for W
SHEAR_N_DESIGN_CASE = 2^9;

SHEAR_H_DESIGN_MM = [300; 600];
SHEAR_AS_DESIGN_MM = [40; 65];
SHEAR_H2B_DESIGN = 2;
SHEAR_S2D_DESIGN = 2;

SHEAR_FC_DESIGN_MPA = [30; 50];

SHEAR_BAR_TYPE_DESIGN = 1;
SHEAR_SD_DESIGN_MM = [6; 10];
SHEAR_SS_DESIGN_MM = 300;
SHEAR_FS_DESIGN_MPA = [235; 335];

SHEAR_FRP_TYPE_DESIGN = 1;
SHEAR_BETA_DESIGN_DEG = [90; 45];
SHEAR_T_FRP_DESIGN_MM = [0.1; 0.5];
% SHEAR_STRENGTH_TYPE = [1; 2];
SHEAR_STRENGTH_TYPE = 2;
SHEAR_W_FRP_DESIGN_MM = 50;
% SHEAR_S_FRP_DESIGN_MM = 300;
SHEAR_RO_FRP_DESIGN = [0.001;0.003];
SHEAR_E_FRP_DESIGN_MPA = [100e3; 200e3];
SHEAR_STRAIN_FRP_DESIGN = [1/100; 1.5/100];
SHEAR_F_FRP_DESIGN_MPA = SHEAR_E_FRP_DESIGN_MPA .* SHEAR_STRAIN_FRP_DESIGN;

SHEAR_FRP_CONFIG_DESIGN = 1;

%% basic information for the establishment of design cases (flexure)
FLEXURE_N_DESIGN_CASE = 2^9;
% geometrical properties
FLEXURE_H_DESIGN_MM = [300; 600];
FLEXURE_AS_DESIGN_MM = [40; 65];
FLEXURE_D_CMP_DESIGN_MM = 35;
FLEXURE_H2B_DESIGN = 2;
FLEXURE_SPAN2H_DESIGN = [8; 12];
FLEXURE_SHEAR2SPAN_DESIGN = 0.5;
FLEXURE_FRP_END_DESIGN_MM = 0;
% concrete properties
FLEXURE_FC_DESIGN_MPA = [30; 50];
% steel properties
FLEXURE_FS_DESIGN_MPA = 335;
FLEXURE_FS_CMP_DESIGN_MPA = 235;
FLEXURE_ES_DESIGN_MPA = 200e3;
FLEXURE_ES_CMP_DESIGN_MPA = 200e3;
FLEXURE_RO_STEEL_DESIGN = [0.5/100; 1/100];
FLEXURE_RO_STEEL_CMP_DESIGN = [0; 0.5/100];
% FRP properties
FLEXURE_FRP_TYPE_DESIGN = 1;
FLEXURE_FRP_CONFIG_DESIGN = 1;
FLEXURE_RO_B_FRP_DESIGN = [0.5;1];
FLEXURE_T_FRP_DESIGN_MM = [0.4; 1.5];
% FLEXURE_T_FRP_DESIGN_MM = [0.2; 0.5];
FLEXURE_E_FRP_DESIGN_MPA = [100e3; 200e3];
FLEXURE_STRAIN_FRP_DESIGN = [1.5/100; 2/100];
% FLEXURE_STRAIN_FRP_DESIGN = [1/100; 1.5/100];

%% Information of random variables
H_NORM_MEAN = 3.05; H_STD_MM = 6.35;

BETA_BIAS = 1; BETA_STD_DEG = 1;

F_FRP_BIAS = 1.00;
% F_FRP_BIAS = 1.42;
F_FRP_COV = 0.12;

E_FRP_BIAS = 1.0;
E_FRP_COV = 0.10;

FC_BIAS = 1.25;
% FC_BIAS = 1/(1-1.645*0.2);
FC_COV = 0.20;
FCT_BIAS = 1.25;
% FCT_BIAS = 1/(1-1.645*0.2);
FCT_COV = 0.20;

FS_BIAS = 1.20;
FS_COV = 0.10;

AS_BIAS = 1.0;
AS_COV = 2.4/100;

LIVE_BIAS = 1.0; LIVE_COV = 0.25;
DEAD_BIAS = 1.05; DEAD_COV = 0.10;

%% Information for reliability analysis
% FACTOR_FRP = (0.50:0.05:1.00)';
FACTOR_FRP = (1.00:0.05:1.50)';
LOAD_RATIO = [0.50;2.0];
ALPHA_MODEL_ERROR = 0.05;
N_MC = 10000;
TARGET_INDEX = 3.0;