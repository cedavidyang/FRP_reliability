% preprocessing based on the user input

%% extract data from .xls(x) files
% B_TEST_ARRAY_MM: beam width (web width for T beams) in mm
% H_TEST_ARRAY_MM: beam height in mm
% D_TEST_ARRAY_MM: beam depth in mm
% DFRP_TEST_ARRAY_MM: FRP depth in mm
% DFRP_TOP_TEST_ARRAY_MM: distance from top of beam to FRP upper-end in mm
% S2D_TEST_ARRAY: shear span ratio
% FC_TEST_ARRAY_MPA: MEAN concrete strength in MPa
% BAR_TYPE_TEST_ARRAY: bar type, 1 for deformed bar
%                                2 for round bar
%                                0 for no bar/not available
% SD_TEST_ARRAY_MM: bar diameter in mm
% SS_TEST_ARRAY_MM: bar interval in mm
% FS_TEST_ARRAY_MPA: bar yielding strength in MPa
% FRP_TYPE_TEST_ARRAY: FRP type, 1 for CFRP
%                                2 for GFRP
%                                0 for others
% FRP_CONFIG_TEST_ARRAY: FRP configuration, 1 for Sheet
%                                           2 for Strip
%                                           0 for error
% FRP_FORM_TEST_ARRAY: FRP strengthening form, 1 for Side bonding
%                                              2 for U-jacketing
%                                              3 for Wrapping
%                                              0 for error
% BETA_TEST_ARRAY_DEG: FRP angle in degree
% T_FRP_TEST_ARRAY_MM: FRP thickness in mm
% E_FRP_TEST_ARRAY_MPA: FRP elastic modulus in MPa
% F_FRP_TEST_ARRAY_MPA: FRP strength in MPa
% W_FRP_TEST_ARRAY_MM: FRP width in mm
% S_FRP_TEST_ARRAY_MM: FRP interval in mm
% V_TOTAL_TEST_ARRAY_KN: test results of total shear resistance in kN
% V_FRP_TEST_ARRAY_KN: FRP contribution to total shear resistance in kN

switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
%         [TEST_DATABASE_NUM, tmpTXT, ~] = xlsread(XLS_DATAFILE_PATH_AND_NAME, 'all');
        load shear_database
        [nRowsXlsNum, ~] = size(TEST_DATABASE_NUM);
        START_ROW = 4; % data starts from START_ROW
        START_COL = 4; % data starts from START_COL
        isNanRows = isnan(TEST_DATABASE_NUM(START_ROW:nRowsXlsNum, START_COL));
        END_ROW = find(isNanRows, 1)+START_ROW-2; % data ends at END_ROW
        END_ROW( isempty(END_ROW) ) = nRowsXlsNum;
        
        
        % geometrical properties
        B_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,5);
        H_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,6);
        D_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,7);
        DFRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,8);
        DFRP_TOP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,9);
        S2D_TEST_ARRAY = TEST_DATABASE_NUM(START_ROW:END_ROW, 12);
        BEAM_SECTION_ARRAY = tmpTXT(START_ROW:END_ROW, 10);
        
        
        % concrete properties
        FC_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW,4);
        
        
        
        % steel properties 
        BAR_TYPE_TEST_ARRAY = zeros(END_ROW - START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 13);
        tmp1 = strcmpi( tmp0, 'D' );
        tmp2 = strcmpi( tmp0, 'R' );
        tmp3 = (~strcmpi(tmp0,'D')) & (~strcmpi(tmp0,'R'));
        BAR_TYPE_TEST_ARRAY( tmp1 ) = 1;
        BAR_TYPE_TEST_ARRAY( tmp2 ) = 2;
        BAR_TYPE_TEST_ARRAY( tmp3 ) = 0;
        clear tmp0 tmp1 tmp2 tmp3
        
        SD_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,14);
        SS_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,15);
        FS_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW,17);
        SD_TEST_ARRAY_MM( BAR_TYPE_TEST_ARRAY==0 ) = 0;
        SS_TEST_ARRAY_MM( BAR_TYPE_TEST_ARRAY==0 ) = 0;
        FS_TEST_ARRAY_MPA( BAR_TYPE_TEST_ARRAY==0 ) = 0;
        
        
        
        % FRP properties
        FRP_TYPE_TEST_ARRAY = zeros(END_ROW - START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 19);
        tmp1 = strncmpi( tmp0, 'CFRP', 4 );
        tmp2 = strncmpi( tmp0, 'GFRP', 4 );
        tmp3 = (~strncmpi(tmp0,'CFRP', 4)) & (~strncmpi(tmp0,'GFRP', 4));
        FRP_TYPE_TEST_ARRAY( tmp1 ) = 1;
        FRP_TYPE_TEST_ARRAY( tmp2 ) = 2;
        FRP_TYPE_TEST_ARRAY( tmp3 ) = 0;
        clear tmp0 tmp1 tmp2 tmp3
        
        FRP_CONFIG_TEST_ARRAY = zeros(END_ROW - START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 20);
        tmp1 = strncmpi( tmp0, 'Sheet', 5 );
        tmp2 = strncmpi( tmp0, 'Strip', 5 );
        tmp3 = (~strncmpi(tmp0,'Sheet', 5)) & (~strncmpi(tmp0,'Strip', 5));
        if sum(tmp3) ~= 0
            disp('[preprocessing]: unknown FRP configuration')
        end
        FRP_CONFIG_TEST_ARRAY( tmp1 ) = 1;
        FRP_CONFIG_TEST_ARRAY( tmp2 ) = 2;
        clear tmp0 tmp1 tmp2 tmp3
        
        FRP_FORM_TEST_ARRAY = zeros(END_ROW - START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 21);
        tmp1 = strcmpi( tmp0, 'Side' );
        tmp2 = strcmpi( tmp0, 'U' );
        tmp3 = strcmpi( tmp0, 'W' );
        tmp4 = (~strcmpi(tmp0,'Side')) & (~strcmpi(tmp0,'U')) & ...
               (~strcmpi(tmp0,'W'));
        if sum(tmp4) ~= 0
            disp('[preprocessing]: unknown strengthening form')
        end
        FRP_FORM_TEST_ARRAY( tmp1 ) = 1;
        FRP_FORM_TEST_ARRAY( tmp2 ) = 2;
        FRP_FORM_TEST_ARRAY( tmp3 ) = 3;
        clear tmp0 tmp1 tmp2 tmp3 tmp4
        
        BETA_TEST_ARRAY_DEG = TEST_DATABASE_NUM(START_ROW:END_ROW,28);
        T_FRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,23);
        E_FRP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW,22);
        F_FRP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW,24);
        W_FRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,25);
        S_FRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,26);
        
        
        
        % test results
        V_TOTAL_TEST_ARRAY_KN = TEST_DATABASE_NUM(START_ROW:END_ROW,30);
        V_FRP_TEST_ARRAY_KN = TEST_DATABASE_NUM(START_ROW:END_ROW,31);
        
        clear tmpTXT
    
    case{'flexure+all', 'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
%         [TEST_DATABASE_NUM, tmpTXT, ~] = xlsread(XLS_DATAFILE_PATH_AND_NAME, 'all');
        load flexure_database
        [nRowsXlsNum, ~] = size(TEST_DATABASE_NUM);
        START_ROW = 4; % data starts from START_ROW
        START_COL = 4; % data starts from START_COL
        isNanRows = isnan(TEST_DATABASE_NUM(START_ROW:nRowsXlsNum, START_COL));
        END_ROW = find(isNanRows, 1)+START_ROW-2; % data ends at END_ROW
        END_ROW( isempty(END_ROW) ) = nRowsXlsNum;
        
        % geometric properties
        B_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 4);
        H_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 5);
        BF_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 6);
        TF_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 7);
        FRP_END_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 8);
        SHEAR_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 9);
        SPAN_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 10);
        D_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 11);
        D_CMP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW,12);
        
        % concrete properties
        FC_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 13);
        
        % steel properties
        ES_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 14)*1e3;
        FS_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 15);
        AREA_STEEL_TEST_ARRAY_MM2 = TEST_DATABASE_NUM(START_ROW:END_ROW, 16);
        ES_CMP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 17)*1e3;
        FS_CMP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 18);
        AREA_STEEL_CMP_TEST_ARRAY_MM2 = TEST_DATABASE_NUM(START_ROW:END_ROW, 19);
        
        % FRP properties
        FRP_TYPE_TEST_ARRAY = zeros(END_ROW-START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 21);
        tmp1 = strncmpi(tmp0, 'C', 1) & (~strncmpi(tmp0, 'CG', 2)) &...
               (~strncmpi(tmp0, 'CA', 2));
        tmp2 = strncmpi(tmp0, 'G', 1);
        FRP_TYPE_TEST_ARRAY(tmp1) = 1;
        FRP_TYPE_TEST_ARRAY(tmp2) = 2;
        FRP_TYPE_TEST_ARRAY( (~tmp1)&(~tmp2) ) = 0;
        
        FRP_CONFIG_TEST_ARRAY = zeros(END_ROW-START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 21);
        tmp1 = strcmpi(tmp0, 'C-W') | strcmpi(tmp0, 'G-W') |...
               strcmpi(tmp0, 'A-W') | strcmpi(tmp0, 'CG-W');        
        tmp2 = strcmpi(tmp0, 'C-P') | strcmpi(tmp0, 'G-P') |...
               strcmpi(tmp0, 'A-P') | strcmpi(tmp0, 'CG-P');
        FRP_CONFIG_TEST_ARRAY(tmp1) = 1;
        FRP_CONFIG_TEST_ARRAY(tmp2) = 2;
        FRP_CONFIG_TEST_ARRAY((~tmp1)&(~tmp2)) = 3;               
        
        E_FRP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 22)*1e3;
        F_FRP_TEST_ARRAY_MPA = TEST_DATABASE_NUM(START_ROW:END_ROW, 23);
        T_FRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 24);
        B_FRP_TEST_ARRAY_MM = TEST_DATABASE_NUM(START_ROW:END_ROW, 25);
        
        % test results
        M_TOTAL_TEST_ARRAY_KNM = TEST_DATABASE_NUM(START_ROW:END_ROW, 27);
        FAIL_MODE_TEST_ARRAY = zeros(END_ROW-START_ROW+1, 1);
        tmp0 = tmpTXT(START_ROW:END_ROW, 29);
        tmp1 = strcmpi(tmp0, 'IC');
        tmp2 = strcmpi(tmp0, 'rupture');
        FAIL_MODE_TEST_ARRAY(tmp1) = 1;
        FAIL_MODE_TEST_ARRAY(tmp2) = 2;
        % assume all the beams failed by FRP have reliable anchors
        ANCHOR_TEST_ARRAY = zeros(END_ROW-START_ROW+1, 1);
        ANCHOR_TEST_ARRAY(tmp1) = 0;
        ANCHOR_TEST_ARRAY(tmp2) = 1;
          
    otherwise
        disp('[preprocessing]: unknown database')
end

%% Design cases
switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
        N_DESIGN_CASE = SHEAR_N_DESIGN_CASE;
    case {'flexure+IC', 'flexure+ic'}
        N_DESIGN_CASE = FLEXURE_N_DESIGN_CASE;
        FLEXURE_ANCHOR_DESIGN = 0;
    case {'flexure+rup', 'flexure+rupture'}
        N_DESIGN_CASE = FLEXURE_N_DESIGN_CASE;
        FLEXURE_ANCHOR_DESIGN = 1;        
    otherwise
end
        
switch SUB_TEST_DATABASE_NAME
    case {'shear+side'}
        FRP_FORM_DESIGN = 1;
    case {'shear+U'}
        FRP_FORM_DESIGN = 2;
    case {'shear+W'}
        FRP_FORM_DESIGN = 3;
    otherwise
end

%% Reliability analysis
switch DESIGN_CODE
    case {'ACI', 'aci'}
        LOAD_FACTOR = [1.4; 1.7];    %PHI=0.85
%         LOAD_FACTOR = [1.2; 1.6];    %PHI=0.75
    case {'HK', 'hk'}
        LOAD_FACTOR = [1.4; 1.6];
    case {'GB', 'gb'}
        LOAD_FACTOR = [1.2; 1.4];
    case {'TR', 'tr', 'FIB', 'fib'}
        LOAD_FACTOR = [1.35; 1.5];      
    case {'ACInew', 'acinew'}
        LOAD_FACTOR = [1.4; 1.7];    %PHI=0.85
%         LOAD_FACTOR = [1.2; 1.6];    %PHI=0.75        
    otherwise
end
        

