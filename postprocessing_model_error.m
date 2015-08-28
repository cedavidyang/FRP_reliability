% post processing of model error analysis

switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
        %% model error diagram based on beam width
        subplot(2,2,1);
        plot(B_TEST_ARRAY_MM, modelError, '.')
        xlabel('Beam width (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{1} = 'Beam width';
        corrModelErrorTest(1) = corr(modelError, B_TEST_ARRAY_MM);
        
        subplot(2,2,2);
        plot(H_TEST_ARRAY_MM, modelError, '.')
        xlabel('Beam height (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{2} = 'Beam height';
        corrModelErrorTest(2) = corr(modelError, H_TEST_ARRAY_MM);
        
        subplot(2,2,3);
        tmpFrpLength = DFRP_TEST_ARRAY_MM - DFRP_TOP_TEST_ARRAY_MM;
        plot(tmpFrpLength, modelError, '.')
        xlabel('FRP length (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{3} = 'FRP length';
        corrModelErrorTest(3) = corr(modelError, tmpFrpLength);
        
        subplot(2,2,4);
        plot(S2D_TEST_ARRAY, modelError, '.')
        xlabel('Shear span ratio')
        ylabel('Model error')
        nameCorrModelErrorTest{4} = 'Shear span ratio';
        corrModelErrorTest(4) = corr(modelError, S2D_TEST_ARRAY);
        
        %% concrete properties
        figure;
        plot(FC_TEST_ARRAY_MPA, modelError, '.')
        xlabel('Concrete strength (MPa)')
        ylabel('Model error')
        nameCorrModelErrorTest{5} = 'Concrete strength';
        corrModelErrorTest(5) = corr(modelError, S2D_TEST_ARRAY);
        
        %% steel properties
        subplot(2,2,1);
        tmpBarType = BAR_TYPE_TEST_ARRAY;
        tmpBarType( BAR_TYPE_TEST_ARRAY==0 ) = [];
        tmpModelError = modelError;
        tmpModelError( BAR_TYPE_TEST_ARRAY==0 ) = [];
        boxplot(tmpModelError, tmpBarType, 'labels', {'Deformed', 'Round'})
        xlabel('Type of steel bars')
        ylabel('Model error')
        
        subplot(2,2,2);
        tmpSd = SD_TEST_ARRAY_MM;
        tmpSd( SD_TEST_ARRAY_MM==0 ) = [];
        tmpModelError = modelError;
        tmpModelError( SD_TEST_ARRAY_MM==0 ) = [];
        plot(tmpSd, tmpModelError, '.')
        xlabel('Diameter of steel bars (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{6} = 'Bar diameter';
        corrModelErrorTest(6) = corr(tmpModelError, tmpSd);
        
        subplot(2,2,3);
        tmpSs = SS_TEST_ARRAY_MM;
        tmpSs( SS_TEST_ARRAY_MM==0 ) = [];
        tmpModelError = modelError;
        tmpModelError( SS_TEST_ARRAY_MM==0 ) = [];
        plot(tmpSs, tmpModelError, '.')
        xlabel('Interval of steel bars (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{7} = 'Bar interval';
        corrModelErrorTest(7) = corr(tmpModelError, tmpSs);
        
        subplot(2,2,4);
        tmpFs = FS_TEST_ARRAY_MPA;
        tmpFs( FS_TEST_ARRAY_MPA==0 ) = [];
        tmpModelError = modelError;
        tmpModelError( FS_TEST_ARRAY_MPA==0 ) = [];
        plot(tmpFs, tmpModelError, '.')
        xlabel('Yield strength (MPa)')
        ylabel('Model error')
        nameCorrModelErrorTest{8} = 'Steel yielding strength';
        corrModelErrorTest(8) = corr(tmpModelError, tmpFs);
        
        %% FRP properties (informatic properties)
        figure;
        subplot(2,2,1);
        boxplot(modelError, FRP_TYPE_TEST_ARRAY, 'grouporder', {'1','2','0'},...
            'labels', {'CFRP', 'GFRP', 'Others'})
        xlabel('Type of FRP')
        ylabel('Model error')
        
        subplot(2,2,2);
        boxplot(modelError, FRP_CONFIG_TEST_ARRAY, 'labels', {'Sheet', 'Strip'});
        xlabel('Configuration of FRP')
        ylabel('Model error')
        
        subplot(2,2,3);
        boxplot(modelError, FRP_FORM_TEST_ARRAY, 'labels', {'Side', 'U', 'W'});
        xlabel('Form of FRP strengthening')
        ylabel('Model error')
        
        subplot(2,2,4);
        tmpBeta = BETA_TEST_ARRAY_DEG;
        tmp_is45or90 = (BETA_TEST_ARRAY_DEG == 45) | ( BETA_TEST_ARRAY_DEG == 90 );
        tmpBeta = tmpBeta(tmp_is45or90);
        tmpModelError = modelError(tmp_is45or90);
        boxplot(tmpModelError, tmpBeta);
        xlabel('Beta of FRP')
        ylabel('Model error')
        
        %% FRP properties (material)
        figure;
        subplot(2, 2, 1);
        plot(T_FRP_TEST_ARRAY_MM, modelError, '.')
        xlabel('Thickness of FRP (mm)')
        ylabel('Model error')
        nameCorrModelErrorTest{9} = 'FRP thickness';
        corrModelErrorTest(9) = corr(modelError, T_FRP_TEST_ARRAY_MM);
        
        subplot(2, 2, 2);
        plot(E_FRP_TEST_ARRAY_MPA, modelError, '.')
        xlabel('Modulus of FRP (MPa)')
        ylabel('Model error')
        nameCorrModelErrorTest{10} = 'FRP modulus';
        corrModelErrorTest(10) = corr(modelError, E_FRP_TEST_ARRAY_MPA);
        
        subplot(2, 2, 3);
        plot(F_FRP_TEST_ARRAY_MPA, modelError, '.')
        xlabel('Strength of FRP (MPa)')
        ylabel('Model error')
        nameCorrModelErrorTest{11} = 'FRP strength';
        corrModelErrorTest(11) = corr(modelError, F_FRP_TEST_ARRAY_MPA);
        
        %% FRP properties (geometrics)
        figure;
        subplot(2,2,1);
        tmpBetaRAD = BETA_TEST_ARRAY_DEG * pi / 180;
        tmp_isContinuous = ( W_FRP_TEST_ARRAY_MM == 1) | ...
            ( W_FRP_TEST_ARRAY_MM == S_FRP_TEST_ARRAY_MM .* sin(tmpBetaRAD) );
        boxplot(modelError, tmp_isContinuous, 'grouporder', {'1','0'},...
            'labels', {'Continuous', 'Discrete'});
        xlabel('Interval of FRP')
        ylabel('Model error')
        
        subplot(2,2,3);
        tmpWidthFRP = W_FRP_TEST_ARRAY_MM( ~tmp_isContinuous );
        tmpModelError = modelError( ~tmp_isContinuous );
        plot(tmpWidthFRP, tmpModelError, '.')
        xlabel('Width of FRP (mm)');
        ylabel('Model error')
        nameCorrModelErrorTest{12} = 'FRP width';
        corrModelErrorTest(12) = corr(tmpModelError, tmpWidthFRP);
        
        subplot(2,2,4);
        tmpIntervalFRP = S_FRP_TEST_ARRAY_MM( ~tmp_isContinuous );
        tmpModelError = modelError( ~tmp_isContinuous );
        plot(tmpIntervalFRP, tmpModelError, '.')
        xlabel('Interval of FRP (mm)');
        ylabel('Model error')
        nameCorrModelErrorTest{13} = 'FRP interval';
        corrModelErrorTest(13) = corr(tmpModelError, tmpIntervalFRP);
        
        figure; bar(corrModelErrorTest);
        set(gca, 'xticklabel', nameCorrModelErrorTest)
        
        isSide = FRP_FORM_TEST_ARRAY == 1;
        isU = FRP_FORM_TEST_ARRAY == 2;
        isW = FRP_FORM_TEST_ARRAY == 3;
        model_error_analysis(resistanceFromPrediction(isSide), resistanceFromTest(isSide));
        model_error_analysis(resistanceFromPrediction(isU), resistanceFromTest(isU));
        model_error_analysis(resistanceFromPrediction(isW), resistanceFromTest(isW));
    case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}        
        figure;
        isIC = FAIL_MODE_TEST_ARRAY ==1;
        isRupture = FAIL_MODE_TEST_ARRAY ==2;
        boxplot(modelError, FAIL_MODE_TEST_ARRAY, 'labels', {'IC debonding', 'FRP rupture'});
        xlabel('Failure mode')
        ylabel('Model error')
        model_error_analysis(resistanceFromPrediction(isIC), resistanceFromTest(isIC));
        model_error_analysis(resistanceFromPrediction(isRupture), resistanceFromTest(isRupture));
    otherwise
end