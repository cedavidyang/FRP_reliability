% main_shear
% main function for time-invariant reliability analysis
clear; %clc;
rng(1); % repeatability
cov_interval = 0.01;
running_type = input(strcat('type of running: 1 for shear+side, 2 for shear+U,',...
'3 for shear+W, 4 for flexure+IC,\n 5 for flexure+rupture, 6 for 1-3,', ...
'7 for 4-5, 8 for 1-5\n'));
fc_type = input('fc type: 1 for fcdet, 2 for 125 bias, 3 for cov10, 4 for 0.10:0.01:0.20 det dc\n');
design_code = input('design guidelines/code\n', 's');

switch running_type
    case 1
        ischeme_start = 1;
        ischeme_end = 1;
    case 2
        ischeme_start = 2;
        ischeme_end = 2;  
    case 3
        ischeme_start = 3;
        ischeme_end = 3;  
    case 4
        ischeme_start = 4;
        ischeme_end = 4; 
    case 5
        ischeme_start = 5;
        ischeme_end = 5;  
    case 6
        ischeme_start = 1;
        ischeme_end = 3;  
    case 7
        ischeme_start = 4;
        ischeme_end = 5;  
    case 8
        ischeme_start = 1;
        ischeme_end = 6; 
    otherwise
        fprintf('illegal running type');
        break;
end

switch fc_type
    case 1
        fccov_start = 0.20;
        fccov_end = 0.20;
        get_bias = @(x) 1/(1-1.645*x);
        fc_name = @(x) 'fcdet';
    case 2
        fccov_start = 0.20;
        fccov_end = 0.20;
        get_bias = @(x) x*0+1.25; 
        fc_name = @(x) 'fc125';
    case 3
        fccov_start = 0.10;
        fccov_end = 0.10;
        get_bias = @(x) 1/(1-1.645*x);
        fc_name = @(x) 'cov10';
    case 4
        fccov_start = 0.10;
        fccov_end = 0.20;
        get_bias = @(x) 1/(1-1.645*x);
        fc_name = @(x) strcat('cov', num2str(100.*x));
    case 5
        fccov_start = 0.15;
        fccov_end = 0.15;
        get_bias = @(x) 1/(1-1.645*x);
        fc_name = @(x) 'cov15';
    otherwise
        fprintf('illegal fc type')
        break
end

switch lower(design_code)
    case {'aci', 'acinew'}
        psi_f = 1.00;
        FACTOR_FRP = (0.10:0.05:2.00)';
        TARGET_INDEX = 3.5;
    case {'hk', 'tr', 'fib'}
        psi_f = (0.10:0.05:2.00)';
        FACTOR_FRP = 1.00;
        TARGET_INDEX = 3.8;
    case {'gb'}
        psi_f = (0.10:0.05:2.00)';
        FACTOR_FRP = 1.00;
        TARGET_INDEX = 3.7;        
end

for ischeme=ischeme_start:ischeme_end
    if fc_type ==4
        re_data = {};
    end
    for fccov=fccov_start:cov_interval:fccov_end;
        if fc_type ~=4
            re_data = {};
        end
        for ipsi = 1:length(psi_f)
            DESIGN_CODE = lower(design_code);
            % material
            FC_COV = fccov;
            FC_BIAS = get_bias(FC_COV);
            FCT_COV = fccov;
            FCT_BIAS = get_bias(FCT_COV);
            switch ischeme
                case 1
                    user_input;
                    SUB_TEST_DATABASE_NAME ='shear+side';
                    % reliability
                    N_MC = 1e4;
                case 2
                    user_input;
                    SUB_TEST_DATABASE_NAME ='shear+U';
                    % reliability
                    N_MC = 1e4;
                case 3
                    user_input;
                    SUB_TEST_DATABASE_NAME ='shear+W';
                    % reliability
                    N_MC = 1e4;
                case 4
                    user_input;
                    SUB_TEST_DATABASE_NAME ='flexure+IC';
                    % reliability
                    N_MC = 500;
                case 5
                    user_input_rupture;
                    SUB_TEST_DATABASE_NAME ='flexure+rupture';
                    % reliability
                    N_MC = 500;
                otherwise
                    fprintf('illegal running type');
            end
            preprocessing;
%             clearvars -except *TEST* *BIAS* *MEAN* *COV* *STD* cov_interval ...
%                 running_type fc_type design_code ischeme_start ischeme_end ...
%                 fccov_start fccov_end get_bias psi_f FACTOR_FRP TARGET_INDEX ...
%                 DESIGN_CODE FC_COV FC_BIAS FCT_COV FCT_BIAS SUB_TEST_DATABASE_NAME ...
%                 N_MC
            save('tmpdata.mat', '*TEST*', '*BIAS*', '*MEAN*', '*COV*', '*STD*');

            %% Model error analysis
            switch DESIGN_CODE
                case {'ACI' 'aci'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, le_warning, ~] = shear_total_ACI('MODEL_ERROR', 1.00);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                        [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_ACI('MODEL_ERROR', 1.00);
                        resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
                    else
                        disp('[main_function]: unknown resistance');
                    end
                case {'hk', 'HK'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, ~] = shear_total_HK('MODEL_ERROR', 1.00);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                        [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_HK('MODEL_ERROR', 1.00);
                        resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
                    else
                        disp('[main_function]: unknown resistance');
                    end  
                case {'GB', 'gb'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, ro, shearReinforceKN] = shear_total_GB('MODEL_ERROR', 1);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                        [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_GB('MODEL_ERROR', 1.00);
                        resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
                    else
                        disp('[main_function]: unknown resistance');
                    end
                case {'TR', 'tr'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, yield_warning, shearReinforceKN] = shear_total_TR('MODEL_ERROR', 1);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                        [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_TR('MODEL_ERROR', 1.00);
                        resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
                    else
                        disp('[main_function]: unknown resistance');
                    end
                case {'FIB', 'fib'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, yield_warning, shearReinforceKN] = shear_total_fib('MODEL_ERROR', 1);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                        [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_fib('MODEL_ERROR', 1.00);
                        resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
                    else
                        disp('[main_function]: unknown resistance');
                    end                    
                case {'ACInew', 'acinew'}
                    if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                        [resistanceFromPrediction, isOverReinforce, yield_warning, shearReinforceKN] = shear_total_ACInew('MODEL_ERROR', 1);
                        resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
                    else
                        disp('[main_function]: unknown resistance');
                    end          
                otherwise
                    disp('[main_function]: unknown guidelines');
            end

            modelError = resistanceFromTest ./ resistanceFromPrediction;
            paramModelError = lognfit(modelError, ALPHA_MODEL_ERROR);
            if ipsi == 1
                beam_resistance = strtok(SUB_TEST_DATABASE_NAME, '+');
                save(strcat('./data/', 'me_', DESIGN_CODE, '_', beam_resistance, '_fccov', num2str(fccov), '.mat'),...
                    'modelError', 'resistanceFromTest', 'resistanceFromPrediction',...
                    'failModeFromPrediction');
            end
            
            % postprocessing_model_error;
            delete 'tmpdata.mat'; pause(3)

            %% Establish design cases
            switch SUB_TEST_DATABASE_NAME
                case {'shear+side', 'shear+U', 'shear+W'}
                    shear_construct_design_case;
                case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
                    flexure_construct_design_case;
                otherwise
            end
            save('tmpdata.mat', '*DESIGN*', '*BIAS*', '*MEAN*', '*COV*', '*STD*');
            nFactorFrp = length(FACTOR_FRP);
            resistanceDesign = zeros(N_DESIGN_CASE, nFactorFrp);
            isOverReinforce = zeros(N_DESIGN_CASE, nFactorFrp);
            roSteel = zeros(N_DESIGN_CASE, nFactorFrp);
            resistReinforce = zeros(N_DESIGN_CASE, nFactorFrp);
            failMode = zeros(N_DESIGN_CASE, nFactorFrp);
            isSteelYielding = zeros(N_DESIGN_CASE, nFactorFrp);

            switch DESIGN_CODE
                case {'ACI' 'aci'}
                    for iFactorFrp = 1:nFactorFrp
                        if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                            [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_ACI('DESIGN_VALUE',...
                                                               FACTOR_FRP(iFactorFrp));
                        elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
                            [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp), isSteelYielding(:,iFactorFrp)]= flexure_total_ACI('DESIGN_VALUE',...
                                                               FACTOR_FRP(iFactorFrp) );
                        else
                            disp('[main_function]: unknown resistance');
                        end
                    end
                case {'HK' 'hk'}
                    for iFactorFrp = 1:nFactorFrp
                        switch SUB_TEST_DATABASE_NAME
                            case {'shear+side', 'shear+U'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_HK('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);  
                            case {'shear+W'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_HK('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);                   
                            case {'flexure+IC', 'flexure+ic'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_HK('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'flexure+rup', 'flexure+rupture'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_HK('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);
                            otherwise
                                disp('[main_function]: unknown resistance');
                        end
                    end
                case {'GB', 'gb'}
                    for iFactorFrp = 1:nFactorFrp
                        switch SUB_TEST_DATABASE_NAME
                            case {'shear+side', 'shear+U'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_GB('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'shear+W'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_GB('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);                    
                            case {'flexure+IC', 'flexure+ic'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_GB('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'flexure+rup', 'flexure+rupture'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_GB('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);
                            otherwise
                                disp('[main_function]: unknown resistance');
                        end
                    end
                case {'TR', 'tr'}
                    for iFactorFrp = 1:nFactorFrp
                        switch SUB_TEST_DATABASE_NAME
                            case {'shear+side', 'shear+U'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_TR('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'shear+W'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_TR('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);
                            case {'flexure+IC', 'flexure+ic'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_TR('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'flexure+rup', 'flexure+rupture'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_TR('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);
                            otherwise
                                disp('[main_function]: unknown resistance');
                        end
                    end                                  
                case {'FIB', 'fib'}
                    for iFactorFrp = 1:nFactorFrp
                        switch SUB_TEST_DATABASE_NAME
                            case {'shear+side', 'shear+U'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_fib('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'shear+W'}
                                [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                    roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_fib('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);                                  
                            case {'flexure+IC', 'flexure+ic'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_fib('DESIGN_VALUE',...
                                    [1.00; FACTOR_FRP(iFactorFrp); psi_f(ipsi)]);
                            case {'flexure+rup', 'flexure+rupture'}
                                [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                                    isSteelYielding(:,iFactorFrp)] = flexure_total_fib('DESIGN_VALUE',...
                                    [FACTOR_FRP(iFactorFrp); 1.00; psi_f(ipsi)]);
                            otherwise
                                disp('[main_function]: unknown resistance');
                        end
                    end  
                case {'ACInew' 'acinew'}
                    for iFactorFrp = 1:nFactorFrp
                        if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
                            [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                                roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_ACInew('DESIGN_VALUE',...
                                                               FACTOR_FRP(iFactorFrp));
                        else
                            disp('[main_function]: unknown resistance');
                        end
                    end        
                otherwise
                    disp('[main_function]: unknown guidelines');
            end
%             delete 'tmpdata.mat'
            % postprocessing_design_case;

%             refine_design_case;

%             save('tmpdata.mat', '*DESIGN*', '*BIAS*', '*MEAN*', '*COV*', '*STD*');
            %% Time-invariant reliability analysis

            nLoadRatio = length(LOAD_RATIO);
            reliabilityResults = zeros(N_DESIGN_CASE, nFactorFrp, nLoadRatio);
            nReliabilityCase = N_DESIGN_CASE*nFactorFrp*nLoadRatio;

            resistMean = zeros(N_DESIGN_CASE, 1);
            resistStd = zeros(N_DESIGN_CASE, 1);
            resistSmp = zeros(N_MC, N_DESIGN_CASE);

            % summary model error
            switch SUB_TEST_DATABASE_NAME
                case {'shear+side', 'shear+U', 'shear+W'}
                    if strcmpi(DESIGN_CODE, 'aci') %&& strcmpi(SUB_TEST_DATABASE_NAME, 'shear+side')
                        modelErrorMean = mean(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN & le_warning~=1));
                        modelErrorStd = std(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN & le_warning~=1));
                    else
                        modelErrorMean = mean(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN));
                        modelErrorStd = std(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN));
                    end
                case {'flexure+IC', 'flexure+ic'}
                    modelErrorMean = mean(modelError(FAIL_MODE_TEST_ARRAY == 1));
                    modelErrorStd = std(modelError(FAIL_MODE_TEST_ARRAY == 1));
                case {'flexure+rup', 'flexure+rupture'}
                    modelErrorMean = mean(modelError(FAIL_MODE_TEST_ARRAY == 2 & failModeFromPrediction==2));
                    modelErrorStd = std(modelError(FAIL_MODE_TEST_ARRAY == 2 & failModeFromPrediction==2)); 
                otherwise
            end
            clearvars -except N_DESIGN_CASE N_MC LOAD_RATIO LOAD_FACTOR LIVE_BIAS LIVE_COV ...
                DEAD_BIAS DEAD_COV SUB_TEST_DATABASE_NAME nFactorFrp nLoadRatio reliabilityResults ...
                nReliabilityCase resistMean resistStd resistSmp modelErrorMean modelErrorStd ...
                resistanceDesign TARGET_INDEX DESIGN_CODE SUB_TEST_DATABASE_NAME ...
                psi_f ipsi ischeme fccov fc_name ...
                cov_interval ...
                running_type fc_type design_code ischeme_start ischeme_end ...
                fccov_start fccov_end get_bias FACTOR_FRP TARGET_INDEX ...
                N_MC re_data
            
            matlabpool 6
            parfor iDesignCase = 1:N_DESIGN_CASE
%             for iDesignCase = 1:N_DESIGN_CASE    
                switch SUB_TEST_DATABASE_NAME
                    case {'shear+side', 'shear+U', 'shear+W'}
                        [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_shear(iDesignCase, N_MC);
                    case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
%                         [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_flexure_rosenblueth(iDesignCase, N_MC);
                        [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_flexure_lhs(iDesignCase, N_MC);
                    otherwise
                end
                resistMean(iDesignCase) = tmpResistMean;
                resistStd(iDesignCase) = tmpResistStd;
                resistSmp(:, iDesignCase) = tmpResistSmp;
                tmpBetaArray = zeros(nFactorFrp, nLoadRatio);
                for jFactorFrp = 1:nFactorFrp
                    for kLoadRatio = 1:nLoadRatio
                        % determine load
                        s = resistanceDesign(iDesignCase, jFactorFrp);
                        ld = LOAD_RATIO(kLoadRatio);
                        liveNorm = ld*s/(LOAD_FACTOR(1)+LOAD_FACTOR(2)*ld);
                        liveMean = liveNorm * LIVE_BIAS;
                        liveStd = liveMean * LIVE_COV;
                        deadNorm = s/(LOAD_FACTOR(1)+LOAD_FACTOR(2)*ld);
                        deadMean = deadNorm * DEAD_BIAS;
                        deadStd = deadMean * DEAD_COV;
                        % compute reliablity indices
            %             reliabilityResults(iDesignCase, jFactorFrp, kLoadRatio) =...
            %                 form_mc(modelErrorMean, modelErrorStd,...
            %                         resistMean, resistStd,...
            %                         liveMean, liveStd, ...
            %                         deadMean, deadStd);
            %             tmpVar = 0.14^2-(tmpResistStd/tmpResistMean)^2;
            %             tmpVar(tmpVar<0) = 1e-8;
            %             modelErrorMean = 1.14; modelErrorStd = modelErrorMean*sqrt(tmpVar);
                        tmpBetaArray(jFactorFrp, kLoadRatio) =...
                            form_mc(modelErrorMean, modelErrorStd,...
                                    tmpResistMean, tmpResistStd,...
                                    liveMean, liveStd, ...
                                    deadMean, deadStd);
                        tmpBeta = tmpBetaArray(jFactorFrp, kLoadRatio);
                    end
                end
                reliabilityResults(iDesignCase, :, :) = tmpBetaArray;
            end
            matlabpool close

            delete tmpdata.mat

            % save data
            if fc_type ~=4 || psi_f(ipsi) == 1
                re_data{end+1} = reliabilityResults;
            else
                continue
            end
            
%             norm_RE = zeros(nFactorFrp,1);
%             mean_RE = zeros(nFactorFrp,1);
%             std_RE = zeros(nFactorFrp,1);
%             upper_RE = zeros(nFactorFrp,1);
%             lower_RE = zeros(nFactorFrp,1);
%             for i_factor = 1:nFactorFrp
%                 RE_col = reliabilityResults(:, i_factor, :);
%                 RE_col = RE_col(:);
%                 RE_col( isnan(RE_col) ) = [];
% 
%                 norm_RE(i_factor) = mean((RE_col-TARGET_INDEX).^2);
%                 mean_RE(i_factor) = mean(RE_col);
%                 std_RE(i_factor) = std(RE_col);
%                 upper_RE(i_factor) = max(RE_col);
%                 lower_RE(i_factor) = min(RE_col);
%             end
            % save(strcat('./figures4_3D/', DESIGN_CODE, '_', SUB_TEST_DATABASE_NAME, num2str(psi_f(ipsi)), '.mat'), 'norm_RE', 'mean_RE', 'std_RE','upper_RE', 'lower_RE');            
            % save(strcat('./data/', lower(DESIGN_CODE), '_', SUB_TEST_DATABASE_NAME, '_', fc_name(fccov), '_', num2str(fccov), '_psi', num2str(psi_f(ipsi)), '.mat'));
            % save(strcat('./data/', 'data_', lower(DESIGN_CODE), '_', SUB_TEST_DATABASE_NAME, '_', num2str(fccov), '.mat'), 'fccov', 'mean_RE');

            % postprocessing_reliability
            % delete tmpUsefulVariables.mat
            % delete tmpAllVariables.mat
        end
        if fc_type ~=4 
            save(strcat('./data/', lower(DESIGN_CODE), '_', SUB_TEST_DATABASE_NAME, '_', fc_name(fccov), num2str(fccov), '.mat'));
        end
    end
    if fc_type == 4
        save(strcat('./data/', lower(DESIGN_CODE), '_', SUB_TEST_DATABASE_NAME, '_fccov.mat'));
    end
end
