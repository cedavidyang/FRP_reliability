clear
rng(1)
%% aci
% aci side
load data_aci_shear+side_fcdet
case_num = ceil(512*rand(1,1)); smp_aci_side = resistSmp(:, case_num);
% aci U
load data_aci_shear+U_fcdet
case_num = ceil(512*rand(1,1)); smp_aci_U = resistSmp(:, case_num);
% aci W
load data_aci_shear+W_fcdet
case_num = ceil(512*rand(1,1)); smp_aci_W = resistSmp(:, case_num);
%% modifided aci
% aci side
load data_acinew_shear+side_fcdet
case_num = ceil(512*rand(1,1)); smp_new_side = resistSmp(:, case_num);
% aci U
load data_acinew_shear+U_fcdet
case_num = ceil(512*rand(1,1)); smp_new_U = resistSmp(:, case_num);
% aci W
load data_acinew_shear+W_fcdet
case_num = ceil(512*rand(1,1)); smp_new_W = resistSmp(:, case_num);