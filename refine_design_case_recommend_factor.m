% refine design cases
switch SUB_TEST_DATABASE_NAME
    case {'flexure+IC', 'flexure+ic'}
        switch DESIGN_CODE
            case {'ACI', 'aci'}
                iCodeFactor = find(FACTOR_FRP==0.85);
            case {'HK', 'hk'}
                iCodeFactor = find(FACTOR_FRP==1.25);
            case {'GB', 'gb'}
                iCodeFactor = find(FACTOR_FRP==1.00);
            otherwise
        end
        tmpFailMode = failMode(:, iCodeFactor);
        fprintf('Number of IC debonding: %d\n', sum(tmpFailMode==1));
        isIC = tmpFailMode==1;
        N_DESIGN_CASE = sum(isIC);
        % geometric properties
        B_DESIGN_ARRAY_MM = B_DESIGN_ARRAY_MM(isIC);
        H_DESIGN_ARRAY_MM = H_DESIGN_ARRAY_MM(isIC);
        BF_DESIGN_ARRAY_MM = BF_DESIGN_ARRAY_MM(isIC);
        TF_DESIGN_ARRAY_MM = TF_DESIGN_ARRAY_MM(isIC);
        FRP_END_DESIGN_ARRAY_MM = FRP_END_DESIGN_ARRAY_MM(isIC);
        SHEAR_DESIGN_ARRAY_MM = SHEAR_DESIGN_ARRAY_MM(isIC);
        SPAN_DESIGN_ARRAY_MM = SPAN_DESIGN_ARRAY_MM(isIC);
        D_DESIGN_ARRAY_MM = D_DESIGN_ARRAY_MM(isIC);
        D_CMP_DESIGN_ARRAY_MM = D_CMP_DESIGN_ARRAY_MM(isIC);
        
        % concrete properties
        FC_DESIGN_ARRAY_MPA = FC_DESIGN_ARRAY_MPA(isIC);
        
        % steel properties
        ES_DESIGN_ARRAY_MPA = ES_DESIGN_ARRAY_MPA(isIC);
        FS_DESIGN_ARRAY_MPA = FS_DESIGN_ARRAY_MPA(isIC);
        AREA_STEEL_DESIGN_ARRAY_MM2 = AREA_STEEL_DESIGN_ARRAY_MM2(isIC);
        ES_CMP_DESIGN_ARRAY_MPA = ES_CMP_DESIGN_ARRAY_MPA(isIC);
        FS_CMP_DESIGN_ARRAY_MPA = FS_CMP_DESIGN_ARRAY_MPA(isIC);
        AREA_STEEL_CMP_DESIGN_ARRAY_MM2 = AREA_STEEL_CMP_DESIGN_ARRAY_MM2(isIC);
        
        % FRP properties
        ANCHOR_DESIGN_ARRAY = ANCHOR_DESIGN_ARRAY(isIC);   
        FRP_TYPE_DESIGN_ARRAY = FRP_TYPE_DESIGN_ARRAY(isIC);
        FRP_CONFIG_DESIGN_ARRAY = FRP_CONFIG_DESIGN_ARRAY(isIC);                 
        E_FRP_DESIGN_ARRAY_MPA = E_FRP_DESIGN_ARRAY_MPA(isIC);
        F_FRP_DESIGN_ARRAY_MPA = F_FRP_DESIGN_ARRAY_MPA(isIC);
        STRAIN_FRP_DESIGN_ARRAY = STRAIN_FRP_DESIGN_ARRAY(isIC);
        T_FRP_DESIGN_ARRAY_MM = T_FRP_DESIGN_ARRAY_MM(isIC);
        B_FRP_DESIGN_ARRAY_MM = B_FRP_DESIGN_ARRAY_MM(isIC);
        
        % test results
        resistanceDesign = resistanceDesign(isIC, :);
        failMode = failMode(isIC, :);
        isSteelYielding = isSteelYielding(isIC, :);  
    case {'flexure+rup', 'flexure+rupture'}
        switch DESIGN_CODE
            case {'ACI', 'aci'}
                iCodeFactor = find(FACTOR_FRP==0.85);
            case {'HK', 'hk'}
                iCodeFactor = find(FACTOR_FRP==1.40);
            case {'GB', 'gb'}
                iCodeFactor = find(FACTOR_FRP==1.40);
            otherwise
        end
        tmpFailMode = failMode(:, iCodeFactor);
        fprintf('Number of FRP rupture: %d\n', sum(tmpFailMode==2));
        isRupture = tmpFailMode==2;
        N_DESIGN_CASE = sum(isRupture);
        % geometric properties
        B_DESIGN_ARRAY_MM = B_DESIGN_ARRAY_MM(isRupture);
        H_DESIGN_ARRAY_MM = H_DESIGN_ARRAY_MM(isRupture);
        BF_DESIGN_ARRAY_MM = BF_DESIGN_ARRAY_MM(isRupture);
        TF_DESIGN_ARRAY_MM = TF_DESIGN_ARRAY_MM(isRupture);
        FRP_END_DESIGN_ARRAY_MM = FRP_END_DESIGN_ARRAY_MM(isRupture);
        SHEAR_DESIGN_ARRAY_MM = SHEAR_DESIGN_ARRAY_MM(isRupture);
        SPAN_DESIGN_ARRAY_MM = SPAN_DESIGN_ARRAY_MM(isRupture);
        D_DESIGN_ARRAY_MM = D_DESIGN_ARRAY_MM(isRupture);
        D_CMP_DESIGN_ARRAY_MM = D_CMP_DESIGN_ARRAY_MM(isRupture);
        
        % concrete properties
        FC_DESIGN_ARRAY_MPA = FC_DESIGN_ARRAY_MPA(isRupture);
        
        % steel properties
        ES_DESIGN_ARRAY_MPA = ES_DESIGN_ARRAY_MPA(isRupture);
        FS_DESIGN_ARRAY_MPA = FS_DESIGN_ARRAY_MPA(isRupture);
        AREA_STEEL_DESIGN_ARRAY_MM2 = AREA_STEEL_DESIGN_ARRAY_MM2(isRupture);
        ES_CMP_DESIGN_ARRAY_MPA = ES_CMP_DESIGN_ARRAY_MPA(isRupture);
        FS_CMP_DESIGN_ARRAY_MPA = FS_CMP_DESIGN_ARRAY_MPA(isRupture);
        AREA_STEEL_CMP_DESIGN_ARRAY_MM2 = AREA_STEEL_CMP_DESIGN_ARRAY_MM2(isRupture);
        
        % FRP properties
        ANCHOR_DESIGN_ARRAY = ANCHOR_DESIGN_ARRAY(isRupture);   
        FRP_TYPE_DESIGN_ARRAY = FRP_TYPE_DESIGN_ARRAY(isRupture);
        FRP_CONFIG_DESIGN_ARRAY = FRP_CONFIG_DESIGN_ARRAY(isRupture);                 
        E_FRP_DESIGN_ARRAY_MPA = E_FRP_DESIGN_ARRAY_MPA(isRupture);
        F_FRP_DESIGN_ARRAY_MPA = F_FRP_DESIGN_ARRAY_MPA(isRupture);
        STRAIN_FRP_DESIGN_ARRAY = STRAIN_FRP_DESIGN_ARRAY(isRupture);
        T_FRP_DESIGN_ARRAY_MM = T_FRP_DESIGN_ARRAY_MM(isRupture);
        B_FRP_DESIGN_ARRAY_MM = B_FRP_DESIGN_ARRAY_MM(isRupture);
        
        % test results
        resistanceDesign = resistanceDesign(isRupture, :);
        failMode = failMode(isRupture, :);
        isSteelYielding = isSteelYielding(isRupture, :);         
    case{'shear+side', 'shear+U', 'shear+W'}
        switch DESIGN_CODE
            case {'ACI', 'aci'}
                iCodeFactor = find(FACTOR_FRP==0.85);
            case {'HK', 'hk'}
                switch SUB_TEST_DATABASE_NAME
                    case{'shear+W'}
                        iCodeFactor = find(FACTOR_FRP==1.40);
                    otherwise
                        iCodeFactor = find(FACTOR_FRP==1.25);
                end
            case {'GB', 'gb'}
                switch SUB_TEST_DATABASE_NAME
                    case{'shear+W'}
                        iCodeFactor = find(FACTOR_FRP==1.40);
                    otherwise
                        iCodeFactor = find(FACTOR_FRP==1.00);
                end
            otherwise
        end        
        tmpisOverReinforce = isOverReinforce(:, iCodeFactor);
        fprintf('Number of diagonal compression: %d\n', sum(tmpisOverReinforce));
        N_DESIGN_CASE = sum(~tmpisOverReinforce);
        % geometrical properties
        B_DESIGN_ARRAY_MM = B_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        H_DESIGN_ARRAY_MM = H_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        D_DESIGN_ARRAY_MM = D_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        AS_DESIGN_ARRAY_MM = AS_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        H2B_DESIGN_ARRAY = H2B_DESIGN_ARRAY(~tmpisOverReinforce);
        DFRP_DESIGN_ARRAY_MM = DFRP_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        DFRP_TOP_DESIGN_ARRAY_MM = DFRP_TOP_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        S2D_DESIGN_ARRAY = S2D_DESIGN_ARRAY(~tmpisOverReinforce);
        % concrete properties
        FC_DESIGN_ARRAY_MPA = FC_DESIGN_ARRAY_MPA(~tmpisOverReinforce);
        % steel properties
        SD_DESIGN_ARRAY_MM = SD_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        FS_DESIGN_ARRAY_MPA = FS_DESIGN_ARRAY_MPA(~tmpisOverReinforce);
        BAR_TYPE_DESIGN_ARRAY = BAR_TYPE_DESIGN_ARRAY(~tmpisOverReinforce);
        SS_DESIGN_ARRAY_MM = SS_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        % FRP properties
        BETA_DESIGN_ARRAY_DEG = BETA_DESIGN_ARRAY_DEG(~tmpisOverReinforce);
        T_FRP_DESIGN_ARRAY_MM = T_FRP_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        E_FRP_DESIGN_ARRAY_MPA = E_FRP_DESIGN_ARRAY_MPA(~tmpisOverReinforce);
        F_FRP_DESIGN_ARRAY_MPA = F_FRP_DESIGN_ARRAY_MPA(~tmpisOverReinforce);
        W_FRP_DESIGN_ARRAY_MM = W_FRP_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        S_FRP_DESIGN_ARRAY_MM = S_FRP_DESIGN_ARRAY_MM(~tmpisOverReinforce);
        STRENGTH_TYPE_ARRAY = STRENGTH_TYPE_ARRAY(~tmpisOverReinforce);
        FRP_TYPE_DESIGN_ARRAY = FRP_TYPE_DESIGN_ARRAY(~tmpisOverReinforce);
        FRP_CONFIG_DESIGN_ARRAY = FRP_CONFIG_DESIGN_ARRAY(~tmpisOverReinforce);
        FRP_FORM_DESIGN_ARRAY = FRP_FORM_DESIGN_ARRAY(~tmpisOverReinforce);  
        % prediction results
        resistanceDesign = resistanceDesign(~tmpisOverReinforce, :);
        isOverReinforce = isOverReinforce(~tmpisOverReinforce, :);
        roSteel = roSteel(~tmpisOverReinforce, :);
        resistReinforce = resistReinforce(~tmpisOverReinforce, :);
    otherwise
end