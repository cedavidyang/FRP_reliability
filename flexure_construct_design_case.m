% construct design cases
iDesignCase = 1;

%% construct design cases
% geometrical properties
H_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
B_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
D_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
SPAN_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
% concrete properties
FC_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
% steel properties
% RO_STEEL_DESIGN_ARRAY = zeros(N_DESIGN_CASE, 1);
% RO_STEEL_CMP_DESIGN_ARRAY = zeros(N_DESIGN_CASE, 1);
AREA_STEEL_DESIGN_ARRAY_MM2 = zeros(N_DESIGN_CASE, 1);
AREA_STEEL_CMP_DESIGN_ARRAY_MM2 = zeros(N_DESIGN_CASE, 1);
% FRP properties
E_FRP_DESIGN_ARRAY_MPA = zeros(N_DESIGN_CASE, 1);
STRAIN_FRP_DESIGN_ARRAY = zeros(N_DESIGN_CASE, 1);
T_FRP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
B_FRP_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);

for iHBeam = 1:length(FLEXURE_H_DESIGN_MM)
    for iSpan = 1:length(FLEXURE_SPAN2H_DESIGN)
        for iFc = 1:length(FLEXURE_FC_DESIGN_MPA)
            for iRoTension = 1:length(FLEXURE_RO_STEEL_DESIGN)
                for iRoCmp = 1:length(FLEXURE_RO_STEEL_CMP_DESIGN)
                    for iRoBFrp = 1:length(FLEXURE_RO_B_FRP_DESIGN)
                        for iTFrp = 1:length(FLEXURE_T_FRP_DESIGN_MM)
                            for iEFrp = 1:length(FLEXURE_E_FRP_DESIGN_MPA)
                                for iStrainFrp = 1:length(FLEXURE_STRAIN_FRP_DESIGN)
                                    % geometrical properties
                                    H_DESIGN_ARRAY_MM(iDesignCase) = FLEXURE_H_DESIGN_MM(iHBeam);
                                    B_DESIGN_ARRAY_MM(iDesignCase) = FLEXURE_H_DESIGN_MM(iHBeam)/FLEXURE_H2B_DESIGN;
                                    D_DESIGN_ARRAY_MM(iDesignCase) = H_DESIGN_ARRAY_MM(iDesignCase) - FLEXURE_AS_DESIGN_MM(iHBeam);
                                    SPAN_DESIGN_ARRAY_MM(iDesignCase) = FLEXURE_SPAN2H_DESIGN(iSpan)*H_DESIGN_ARRAY_MM(iDesignCase) ;
                                    % concrete properties
                                    FC_DESIGN_ARRAY_MPA(iDesignCase) = FLEXURE_FC_DESIGN_MPA(iFc);
                                    % steel properties
                                    AREA_STEEL_DESIGN_ARRAY_MM2(iDesignCase) = B_DESIGN_ARRAY_MM(iDesignCase) * D_DESIGN_ARRAY_MM(iDesignCase) *...
                                                                  FLEXURE_RO_STEEL_DESIGN(iRoTension);
                                    AREA_STEEL_CMP_DESIGN_ARRAY_MM2(iDesignCase) = B_DESIGN_ARRAY_MM(iDesignCase) * D_DESIGN_ARRAY_MM(iDesignCase) *...
                                                                  FLEXURE_RO_STEEL_CMP_DESIGN(iRoCmp);
                                    % FRP properties
                                    E_FRP_DESIGN_ARRAY_MPA(iDesignCase) = FLEXURE_E_FRP_DESIGN_MPA(iEFrp);
                                    STRAIN_FRP_DESIGN_ARRAY(iDesignCase) = FLEXURE_STRAIN_FRP_DESIGN(iStrainFrp);
                                    T_FRP_DESIGN_ARRAY_MM(iDesignCase) = FLEXURE_T_FRP_DESIGN_MM(iTFrp);
                                    B_FRP_DESIGN_ARRAY_MM(iDesignCase) = FLEXURE_RO_B_FRP_DESIGN(iRoBFrp) * B_DESIGN_ARRAY_MM(iDesignCase);
                                    
                                    iDesignCase = iDesignCase + 1;
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
BF_DESIGN_ARRAY_MM = B_DESIGN_ARRAY_MM;
TF_DESIGN_ARRAY_MM = zeros(N_DESIGN_CASE, 1);
FRP_END_DESIGN_ARRAY_MM = FLEXURE_FRP_END_DESIGN_MM * ones(N_DESIGN_CASE, 1);
SHEAR_DESIGN_ARRAY_MM = SPAN_DESIGN_ARRAY_MM * FLEXURE_SHEAR2SPAN_DESIGN;
D_CMP_DESIGN_ARRAY_MM = FLEXURE_D_CMP_DESIGN_MM * ones(N_DESIGN_CASE, 1);
% steel properties
ES_DESIGN_ARRAY_MPA = FLEXURE_ES_DESIGN_MPA * ones(N_DESIGN_CASE, 1);
FS_DESIGN_ARRAY_MPA = FLEXURE_FS_DESIGN_MPA * ones(N_DESIGN_CASE, 1);
ES_CMP_DESIGN_ARRAY_MPA = FLEXURE_ES_CMP_DESIGN_MPA * ones(N_DESIGN_CASE, 1);
FS_CMP_DESIGN_ARRAY_MPA = FLEXURE_FS_CMP_DESIGN_MPA * ones(N_DESIGN_CASE, 1);
% FRP properties
FRP_TYPE_DESIGN_ARRAY = FLEXURE_FRP_TYPE_DESIGN * ones(N_DESIGN_CASE, 1);
FRP_CONFIG_DESIGN_ARRAY = FLEXURE_FRP_CONFIG_DESIGN * ones(N_DESIGN_CASE, 1);
F_FRP_DESIGN_ARRAY_MPA = E_FRP_DESIGN_ARRAY_MPA .* STRAIN_FRP_DESIGN_ARRAY;
ANCHOR_DESIGN_ARRAY = FLEXURE_ANCHOR_DESIGN * ones(N_DESIGN_CASE, 1);