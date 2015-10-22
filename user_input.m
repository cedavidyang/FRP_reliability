% user input file
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

F_FRP_COV = 0.12;
% F_FRP_COV = 0.06;
% F_FRP_BIAS = 1.42;
F_FRP_BIAS = 1/(1-1.645*F_FRP_COV);

E_FRP_BIAS = 1.0;
E_FRP_COV = 0.10;
% E_FRP_COV = 0.05;

FS_COV = 0.10;
% FS_COV = 0.05;
FS_BIAS = 1/(1-1.645*FS_COV);

AS_BIAS = 1.0;
AS_COV = 2.4/100;

LIVE_BIAS = 1.0; LIVE_COV = 0.25;
DEAD_BIAS = 1.05; DEAD_COV = 0.10;

%% Information for reliability analysis
LOAD_RATIO = [0.50;2.0];
ALPHA_MODEL_ERROR = 0.05;