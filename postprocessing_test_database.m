% post processing of test database

%% geometrical properties
switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
        subplot(2,2,1);
        hist(B_TEST_ARRAY_MM);
        xlabel('Beam width (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        hist(H_TEST_ARRAY_MM);
        xlabel('Beam height (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmpFrpLength = DFRP_TEST_ARRAY_MM - DFRP_TOP_TEST_ARRAY_MM;
        hist(tmpFrpLength);
        xlabel('FRP height (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,4);
        hist(S2D_TEST_ARRAY);
        xlabel('shear span ratio')
        ylabel('Number in test database')
        
        %% concrete properties
        figure;
        hist(FC_TEST_ARRAY_MPA);
        xlabel('Concrete strength (MPa)')
        ylabel('Number in test database')
        
        %% steel properties
        figure;
        subplot(2,2,1);
        tmp_nDeformed = sum( BAR_TYPE_TEST_ARRAY == 1);
        tmp_nRound = sum( BAR_TYPE_TEST_ARRAY == 2);
        bar([tmp_nDeformed, tmp_nRound]);
        set(gca, 'xticklabel', {'Deformed', 'Round'});
        xlabel('Type of steel bars')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        tmpSd = SD_TEST_ARRAY_MM;
        tmpSd( SD_TEST_ARRAY_MM==0 ) = [];
        hist(tmpSd);
        xlabel('Diameter of steel bars (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmpSs = SS_TEST_ARRAY_MM;
        tmpSs( SS_TEST_ARRAY_MM==0 ) = [];
        tmpDepthBeam = D_TEST_ARRAY_MM;
        tmpDepthBeam( SS_TEST_ARRAY_MM==0 ) = [];
        hist(tmpSs./tmpDepthBeam);
        xlabel('Ratio of stirrup interval to beam depth')
        ylabel('Number in test database')
        
        subplot(2,2,4);
        tmpFs = FS_TEST_ARRAY_MPA;
        tmpFs( FS_TEST_ARRAY_MPA==0 ) = [];
        hist(tmpFs);
        xlabel('Yield strength (MPa)')
        ylabel('Number in test database')
        
        %% FRP properties
        figure;
        subplot(2,2,1);
        tmp_nCFRP = sum( FRP_TYPE_TEST_ARRAY == 1);
        tmp_nGFRP = sum( FRP_TYPE_TEST_ARRAY == 2);
        tmp_nOtherFRP = sum( FRP_TYPE_TEST_ARRAY == 0);
        bar([tmp_nCFRP, tmp_nGFRP, tmp_nOtherFRP]);
        set(gca, 'xticklabel', {'CFRP', 'GFRP', 'Others'});
        xlabel('Type of FRP')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        tmp_nSheet = sum( FRP_CONFIG_TEST_ARRAY == 1);
        tmp_nStrip = sum( FRP_CONFIG_TEST_ARRAY == 2);
        bar([tmp_nSheet, tmp_nStrip]);
        set(gca, 'xticklabel', {'Sheet', 'Strip'});
        xlabel('Configuration of FRP')
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmp_nSide = sum( FRP_FORM_TEST_ARRAY == 1);
        tmp_nU = sum( FRP_FORM_TEST_ARRAY == 2);
        tmp_nW = sum( FRP_FORM_TEST_ARRAY == 3);
        bar([tmp_nSide, tmp_nU, tmp_nW]);
        set(gca, 'xticklabel', {'Side', 'U', 'W'});
        xlabel('Form of FRP strengthening')
        ylabel('Number in test database')
        
        subplot(2,2,4);
        tmp_nBeta45 = sum( BETA_TEST_ARRAY_DEG == 45);
        tmp_nBeta90 = sum( BETA_TEST_ARRAY_DEG == 90);
        bar([tmp_nBeta45, tmp_nBeta90]);
        set(gca, 'xticklabel', {'45', '90'});
        xlabel('Beta of FRP')
        ylabel('Number in test database')
        
        %% FRP properties (material)
        figure;
        subplot(2, 2, 1);
        hist(T_FRP_TEST_ARRAY_MM)
        xlabel('Thickness of FRP (mm)')
        ylabel('Number in test database')
        
        subplot(2, 2, 2);
        hist(E_FRP_TEST_ARRAY_MPA)
        xlabel('Modulus of FRP (MPa)')
        ylabel('Number in test database')
        
        subplot(2, 2, 3);
        hist(F_FRP_TEST_ARRAY_MPA)
        xlabel('Strength of FRP (MPa)')
        ylabel('Number in test database')
        
        %% FRP properties (geometrics)
        figure;
        subplot(2,2,1);
        tmpBetaRAD = BETA_TEST_ARRAY_DEG * pi / 180;
        tmp_isContinuous = ( W_FRP_TEST_ARRAY_MM == 1) | ...
            ( W_FRP_TEST_ARRAY_MM == S_FRP_TEST_ARRAY_MM .* sin(tmpBetaRAD) );
        tmp_nContinuous = sum(tmp_isContinuous);
        tmp_nDiscrete = sum( ~tmp_isContinuous);
        bar([tmp_nContinuous, tmp_nDiscrete]);
        set(gca, 'xticklabel', {'Continous', 'Discrete'});
        xlabel('Interval of FRP')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        tmpWidthFRP = W_FRP_TEST_ARRAY_MM( ~tmp_isContinuous );
        hist(tmpWidthFRP, 8)
        xlabel('Width of FRP (mm)');
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmpIntervalFRP = S_FRP_TEST_ARRAY_MM( ~tmp_isContinuous );
        tmpNetIntervalFrp = tmpIntervalFRP - tmpWidthFRP;
        tmpBeamDepth = D_TEST_ARRAY_MM( ~tmp_isContinuous );
        hist(tmpIntervalFRP./tmpBeamDepth)
        xlabel('Ratio of FRP interval to beam depth');
        ylabel('Number in test database')
        
        subplot(2,2,4);
        tmpIntervalFRP = S_FRP_TEST_ARRAY_MM( ~tmp_isContinuous );
        tmpRoFrp = 2*T_FRP_TEST_ARRAY_MM(~tmp_isContinuous ).*tmpWidthFRP ./ ...
            B_TEST_ARRAY_MM( ~tmp_isContinuous )./ tmpIntervalFRP;
        hist(tmpRoFrp)
        xlabel('Rinforcement ratio of FRP (mm)');
        ylabel('Number in test database')
    case{'flexure+all', 'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
        %% geometric properties      
        subplot(2,2,1);
        hist(H_TEST_ARRAY_MM);
        xlabel('Beam height (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        hist(SPAN_TEST_ARRAY_MM);
        xlabel('Beam span (mm)')
        ylabel('Number in test database') 
        
        subplot(2,2,3);
        hist(H_TEST_ARRAY_MM./B_TEST_ARRAY_MM);
        xlabel('Aspect ratios of beam sections')
        ylabel('Number in test database') 
        
        subplot(2,2,4);
        hist(SPAN_TEST_ARRAY_MM./H_TEST_ARRAY_MM);
        xlabel('Ratio of span to height')
        ylabel('Number in test database')
                
        %% concrete properties
        figure;
        hist(FC_TEST_ARRAY_MPA);
        xlabel('Concrete strength (MPa)')
        ylabel('Number in test database')
        
        %% steel properties
        figure;
        subplot(2,2,1);
        hist(FS_TEST_ARRAY_MPA);
        xlabel('Strength of tensile steel (MPa)')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        tmpisWithCmp = AREA_STEEL_CMP_TEST_ARRAY_MM2~=0;
        tmpFsCmp = FS_CMP_TEST_ARRAY_MPA(tmpisWithCmp);
        hist(tmpFsCmp);
        xlabel('Strength of compressive steel (MPa)')
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmpRoTensile = AREA_STEEL_TEST_ARRAY_MM2 ./ (B_TEST_ARRAY_MM.*D_TEST_ARRAY_MM);
        hist(tmpRoTensile);
        xlabel('Tensile reinforcement ratio')
        ylabel('Number in test database')
        
        subplot(2,2,4);
        tmpRoCmp = AREA_STEEL_CMP_TEST_ARRAY_MM2(tmpisWithCmp) ./ ...
                   (B_TEST_ARRAY_MM(tmpisWithCmp) .* D_TEST_ARRAY_MM(tmpisWithCmp));
        hist(tmpRoCmp);
        xlabel('Compressive reinforcement ratio')
        ylabel('Number in test database')    

        %% FRP properties (informatics)
        figure;
        subplot(2,2,1);
        tmp_nCFRP = sum( FRP_TYPE_TEST_ARRAY == 1);
        tmp_nGFRP = sum( FRP_TYPE_TEST_ARRAY == 2);
        tmp_nOtherFRP = sum( FRP_TYPE_TEST_ARRAY == 0);
        bar([tmp_nCFRP, tmp_nGFRP, tmp_nOtherFRP]);
        set(gca, 'xticklabel', {'CFRP', 'GFRP', 'Others'});
        xlabel('Type of FRP')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        tmp_nSheet = sum( FRP_CONFIG_TEST_ARRAY == 1);
        tmp_nStrip = sum( FRP_CONFIG_TEST_ARRAY == 2);
        bar([tmp_nSheet, tmp_nStrip]);
        set(gca, 'xticklabel', {'Sheet', 'Strip'});
        xlabel('Configuration of FRP')
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmp_nIc = sum( FAIL_MODE_TEST_ARRAY == 1);
        tmp_nRupture = sum( FAIL_MODE_TEST_ARRAY == 2);
        bar([tmp_nIc, tmp_nRupture]);
        set(gca, 'xticklabel', {'IC debonding', 'FRP rupture'});
        xlabel('Failure mode')
        ylabel('Number in test database')        
        
        %% FRP properties (material)
        figure;        
        subplot(2, 2, 1);
        hist(E_FRP_TEST_ARRAY_MPA)
        xlabel('Modulus of FRP (MPa)')
        ylabel('Number in test database')
        
        subplot(2, 2, 2);
        hist(F_FRP_TEST_ARRAY_MPA)
        xlabel('Strength of FRP (MPa)')
        ylabel('Number in test database')
        
        subplot(2, 2, 3);
        hist(F_FRP_TEST_ARRAY_MPA./E_FRP_TEST_ARRAY_MPA)
        xlabel('Ultimate strain of FRP')
        ylabel('Number in test database')
        
        %% FRP properties (geometrics)
        figure;
        subplot(2, 2, 1);
        hist(T_FRP_TEST_ARRAY_MM)
        xlabel('Thickness of FRP (mm)')
        ylabel('Number in test database')
        
        subplot(2,2,2);
        hist(B_FRP_TEST_ARRAY_MM)
        xlabel('Width of FRP (mm)');
        ylabel('Number in test database')
        
        subplot(2,2,3);
        tmpRoFrp = B_FRP_TEST_ARRAY_MM .* T_FRP_TEST_ARRAY_MM ./...
                   (B_TEST_ARRAY_MM.*H_TEST_ARRAY_MM);
        hist(tmpRoFrp(FAIL_MODE_TEST_ARRAY == 1))
        xlabel('Strengthening ratio of FRP');
        ylabel('Number in test database')
    otherwise
end