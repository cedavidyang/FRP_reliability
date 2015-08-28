% main function for time-invariant reliability analysis

user_input;
% user_input_rupture;
preprocessing;
save('tmpdata.mat', '*TEST*', '*BIAS*', '*MEAN*', '*COV*', '*STD*');

% postprocessing_test_database;

%% Model error analysis

switch DESIGN_CODE
    case {'ACI' 'aci'}
        if strncmpi( SUB_TEST_DATABASE_NAME, 'shear', 5)
            [resistanceFromPrediction, isOverReinforce, ~] = shear_total_ACI('MODEL_ERROR', 1.00);
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
            [resistanceFromPrediction, isOverReinforce, ~] = shear_total_GB('MODEL_ERROR', 1);
            resistanceFromTest = V_TOTAL_TEST_ARRAY_KN;
        elseif strncmpi( SUB_TEST_DATABASE_NAME, 'flexure', 7)
            [resistanceFromPrediction,failModeFromPrediction,~] = flexure_total_GB('MODEL_ERROR', 1.00);
            resistanceFromTest = M_TOTAL_TEST_ARRAY_KNM;
        else
            disp('[main_function]: unknown resistance');
        end        
    otherwise
        disp('[main_function]: unknown guidelines');
end

modelError = resistanceFromTest ./ resistanceFromPrediction;
paramModelError = lognfit(modelError, ALPHA_MODEL_ERROR);

% postprocessing_model_error;
delete 'tmpdata.mat'
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
                case {'shear+side', 'shear+U', 'shear+W'}
                    [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                        roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_HK('DESIGN_VALUE',...
                        [1.40; FACTOR_FRP(iFactorFrp)]);
                case {'flexure+IC', 'flexure+ic'}
                    [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                        isSteelYielding(:,iFactorFrp)] = flexure_total_HK('DESIGN_VALUE',...
                        [1.40; FACTOR_FRP(iFactorFrp)]);
                case {'flexure+rup', 'flexure+rupture'}
                    [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                        isSteelYielding(:,iFactorFrp)] = flexure_total_HK('DESIGN_VALUE',...
                        [FACTOR_FRP(iFactorFrp); 1.25]);
                otherwise
                    disp('[main_function]: unknown resistance');
            end
        end
    case {'GB', 'gb'}
        for iFactorFrp = 1:nFactorFrp
            switch SUB_TEST_DATABASE_NAME
                case {'shear+side', 'shear+U', 'shear+W'}
                    [resistanceDesign(:,iFactorFrp), isOverReinforce(:,iFactorFrp),...
                        roSteel(:,iFactorFrp), resistReinforce(:,iFactorFrp)] = shear_total_GB('DESIGN_VALUE',[FACTOR_FRP(iFactorFrp); 1.25]);
                case {'flexure+IC', 'flexure+ic'}
                    [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                        isSteelYielding(:,iFactorFrp)] = flexure_total_GB('DESIGN_VALUE',...
                        [1.40; FACTOR_FRP(iFactorFrp)]);
                case {'flexure+rup', 'flexure+rupture'}
                    [resistanceDesign(:,iFactorFrp), failMode(:,iFactorFrp),...
                        isSteelYielding(:,iFactorFrp)] = flexure_total_GB('DESIGN_VALUE',...
                        [FACTOR_FRP(iFactorFrp); 1.00]);
                otherwise
                    disp('[main_function]: unknown resistance');
            end
        end
    otherwise
        disp('[main_function]: unknown guidelines');
end
delete 'tmpdata.mat'
% postprocessing_design_case;

refine_design_case;

save('tmpdata.mat', '*DESIGN*', '*BIAS*', '*MEAN*', '*COV*', '*STD*');
%% Time-invariant reliability analysis
save('tmpAllVariables.mat');
nLoadRatio = length(LOAD_RATIO);
reliabilityResults = zeros(N_DESIGN_CASE, nFactorFrp, nLoadRatio);
nReliabilityCase = N_DESIGN_CASE*nFactorFrp*nLoadRatio;

resistMean = zeros(N_DESIGN_CASE, 1);
resistStd = zeros(N_DESIGN_CASE, 1);
resistSmp = zeros(N_MC, N_DESIGN_CASE);

% summary model error
switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
        modelErrorMean = mean(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN));
        modelErrorStd = std(modelError(FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN));
    case {'flexure+IC', 'flexure+ic'}
        modelErrorMean = mean(modelError(FAIL_MODE_TEST_ARRAY == 1));
        modelErrorStd = std(modelError(FAIL_MODE_TEST_ARRAY == 1));
    case {'flexure+rup', 'flexure+rupture'}
        modelErrorMean = mean(modelError(FAIL_MODE_TEST_ARRAY == 2));
        modelErrorStd = std(modelError(FAIL_MODE_TEST_ARRAY == 2)); 
    otherwise
end

save('tmpUsefulVariables.mat', 'N_DESIGN_CASE', 'N_MC', 'LOAD_RATIO', 'LOAD_FACTOR', 'LIVE_BIAS', ...
    'LIVE_COV', 'DEAD_BIAS', 'DEAD_COV', 'SUB_TEST_DATABASE_NAME', ...
    'nFactorFrp', 'nLoadRatio', 'reliabilityResults', 'nReliabilityCase', ...
    'resistMean', 'resistStd', 'resistSmp', 'modelErrorMean', 'modelErrorStd', ...
    'resistanceDesign');
clear
load tmpUsefulVariables

matlabpool 4
parfor iDesignCase = 1:N_DESIGN_CASE
% for iDesignCase = 1:N_DESIGN_CASE    
% for iDesignCase = 458  
    switch SUB_TEST_DATABASE_NAME
        case {'shear+side', 'shear+U', 'shear+W'}
            [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_shear(iDesignCase, N_MC);
        case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
            [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_flexure_rosenblueth(iDesignCase, N_MC);
        otherwise
    end
%     [tmpResistMean, tmpResistStd, tmpResistSmp] = determine_resistance(iDesignCase, N_MC);
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
            % conpute reliablity indices
%             reliabilityResults(iDesignCase, jFactorFrp, kLoadRatio) =...
%                 form_mc(modelErrorMean, modelErrorStd,...
%                         resistMean, resistStd,...
%                         liveMean, liveStd, ...
%                         deadMean, deadStd);
            tmpBetaArray(jFactorFrp, kLoadRatio) =...
                form_mc(modelErrorMean, modelErrorStd,...
                        tmpResistMean, tmpResistStd,...
                        liveMean, liveStd, ...
                        deadMean, deadStd);
            tmpBeta = tmpBetaArray(jFactorFrp, kLoadRatio);
%             if tmpBeta <2.95
%                 fprintf('suspiciously small beta=%f, Rcov = %f, iCase=%d\n',...
%                      [tmpBeta, tmpResistStd/tmpResistMean, iDesignCase]);
%             end
        end
    end
    reliabilityResults(iDesignCase, :, :) = tmpBetaArray;
end
matlabpool close

% resistMean = resistMean';
% resistStd= resistStd';
    
delete tmpdata.mat

load tmpAllVariables
postprocessing_reliability
delete tmpUsefulVariables.mat
delete tmpAllVariables.mat
            