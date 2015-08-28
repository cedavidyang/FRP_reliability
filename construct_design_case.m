% construct design cases
iDesignCase = 1;

%% construct design cases
% geometrical properties
B_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
H_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
D_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
% concrete properties
FC_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
% steel properties 
AS_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
SD_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
FS_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
% FRP properties
BETA_DESIGN_ARRAY_DEG = zeros(N_DESIGN_CASE, 1);
T_FRP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
E_FRP_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
F_FRP_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
W_FRP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
S_FRP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
STRENGTH_TYPE_ARRAY = zeros(N_DESIGN_CASE, 1);

for iFcDesign = 1:length(FC_DESIGN_MPA)
    for iHDesign = 1:length(H_DESIGN_MM)
        for iSdDesign = 1:length(SD_DESIGN_MM)
            for iFsDesign = 1:length(FS_DESIGN_MPA)
                for iBetaDesign = 1:length(BETA_DESIGN_DEG)
                    for iEfDesign = 1:length(E_FRP_DESIGN_MPA)
                        for iStrainDesign = 1:length(STRAIN_FRP_DESIGN)
                            for iTfDesign = 1:length(T_FRP_DESIGN_MM)
%                                 for iStrengthTypeDesign = 1:length(STRENGTH_TYPE)
                                for iRoFRP = 1:length(RO_FRP_DESIGN)
                                    FC_DESIGN_ARRAY_MPA(iDesignCase) = FC_DESIGN_MPA(iFcDesign);
                                    H_DESIGN_ARRAY_MM(iDesignCase) = H_DESIGN_MM(iHDesign);
                                    B_DESIGN_ARRAY_MM(iDesignCase) = H_DESIGN_MM(iHDesign) / H2B_DESIGN;
                                    AS_DESIGN_ARRAY_MM(iDesignCase) = AS_DESIGN_MM(iHDesign);
                                    D_DESIGN_ARRAY_MM(iDesignCase) = H_DESIGN_MM(iHDesign) - AS_DESIGN_MM(iHDesign);
                                    SD_DESIGN_ARRAY_MM(iDesignCase) = SD_DESIGN_MM(iSdDesign);
                                    FS_DESIGN_ARRAY_MPA(iDesignCase) = FS_DESIGN_MPA(iFsDesign);
                                    BETA_DESIGN_ARRAY_DEG(iDesignCase) = BETA_DESIGN_DEG(iBetaDesign);
                                    E_FRP_DESIGN_ARRAY_MPA(iDesignCase) = E_FRP_DESIGN_MPA(iEfDesign);
                                    T_FRP_DESIGN_ARRAY_MM(iDesignCase) = T_FRP_DESIGN_MM(iTfDesign);
                                    F_FRP_DESIGN_ARRAY_MPA(iDesignCase) = E_FRP_DESIGN_MPA(iEfDesign) * STRAIN_FRP_DESIGN(iStrainDesign);
                                    
%                                     W_FRP_DESIGN_ARRAY_MM(iDesignCase) = (STRENGTH_TYPE(iStrengthTypeDesign)==1)*1 +...
%                                         (STRENGTH_TYPE(iStrengthTypeDesign)==2)*W_FRP_DESIGN_MM;
%                                     S_FRP_DESIGN_ARRAY_MM(iDesignCase) = (STRENGTH_TYPE(iStrengthTypeDesign)==1)/sin(BETA_DESIGN_DEG(iBetaDesign)/180*pi) +...
%                                         (STRENGTH_TYPE(iStrengthTypeDesign)==2)*(2*T_FRP_DESIGN_MM(iTfDesign)*W_FRP_DESIGN_MM/...
%                                          B_DESIGN_ARRAY_MM(iDesignCase)/RO_FRP_DESIGN);
                                    roFrpMax = 2*T_FRP_DESIGN_ARRAY_MM(iDesignCase)/B_DESIGN_ARRAY_MM(iDesignCase)*...
                                               sin(BETA_DESIGN_DEG(iBetaDesign)/180*pi);
                                    isContinuous = roFrpMax<RO_FRP_DESIGN(iRoFRP);
                                    STRENGTH_TYPE_ARRAY(iDesignCase) = isContinuous + (~isContinuous)*STRENGTH_TYPE;
                                    W_FRP_DESIGN_ARRAY_MM(iDesignCase) = isContinuous + (~isContinuous)*W_FRP_DESIGN_MM;
                                    roFrp = min(roFrpMax, RO_FRP_DESIGN(iRoFRP));
                                    S_FRP_DESIGN_ARRAY_MM(iDesignCase) = isContinuous/sin(BETA_DESIGN_DEG(iBetaDesign)/180*pi) +...
                                        (~isContinuous) * 2*T_FRP_DESIGN_MM(iTfDesign)*W_FRP_DESIGN_MM/...
                                        B_DESIGN_ARRAY_MM(iDesignCase)/roFrp;
                                    
                                    iDesignCase = iDesignCase+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% other properties
% geometric properties
H2B_DESIGN_ARRAY = H2B_DESIGN*ones(N_DESIGN_CASE, 1);
DFRP_DESIGN_ARRAY_MM = H_DESIGN_ARRAY_MM;
DFRP_TOP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
S2D_DESIGN_ARRAY = S2D_DESIGN * ones(N_DESIGN_CASE,1);
% steel properties
BAR_TYPE_DESIGN_ARRAY = BAR_TYPE_DESIGN * ones(N_DESIGN_CASE,1);
SS_DESIGN_ARRAY_MM = SS_DESIGN_MM * ones(N_DESIGN_CASE,1);
% FRP properties
FRP_TYPE_DESIGN_ARRAY = FRP_TYPE_DESIGN * ones(N_DESIGN_CASE, 1);
FRP_CONFIG_DESIGN_ARRAY = FRP_CONFIG_DESIGN * ones(N_DESIGN_CASE,1);
FRP_FORM_DESIGN_ARRAY = FRP_FORM_DESIGN * ones(N_DESIGN_CASE,1);

% W_FRP_DESIGN_ARRAY_MM = W_FRP_DESIGN_MM * ones(N_DESIGN_CASE,1);
% STRENGTH_TYPE_ARRAY = STRENGTH_TYPE * ones(N_DESIGN_CASE,1);