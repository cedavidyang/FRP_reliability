% post processing of design cases

switch SUB_TEST_DATABASE_NAME
    case {'shear+side', 'shear+U', 'shear+W'}
        iNoFactorFrp = find(FACTOR_FRP==1);
        figure;
        plot(resistanceDesign(:, iNoFactorFrp), 'o');
        hold on
        attention = find(isOverReinforce(:, iNoFactorFrp)==1);
        plot(attention, resistanceDesign(attention, iNoFactorFrp), 'ro');
        figure;
        subplot(1,2,1);
        hist(resistanceDesign(:, iNoFactorFrp));
        xlabel('Total resistance')
        ylabel('Number of design cases');
        subplot(1,2,2);
        switch DESIGN_CODE
            case {'ACI', 'aci'}
                hist( (1-roSteel(:, iNoFactorFrp)).*resistReinforce(:, iNoFactorFrp) ./ resistanceDesign(:, iNoFactorFrp) );
            case {'HK', 'hk'}
                hist( (1-roSteel(:, iNoFactorFrp)).*resistReinforce(:, iNoFactorFrp) ./ resistanceDesign(:, iNoFactorFrp) );
            case {'GB', 'gb'}
                hist( (1-roSteel(:, iNoFactorFrp)).*resistReinforce(:, iNoFactorFrp) ./ resistanceDesign(:, iNoFactorFrp) );
        end
        xlabel('FRP contribution')
        ylabel('Number of design cases');
    case {'flexure+all', 'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
        switch DESIGN_CODE
            case {'ACI', 'aci'}
                iCodeFactor = find(FACTOR_FRP==0.85);
            case {'HK', 'hk'}
            case {'GB', 'gb'}
            otherwise
        end
        tmpFailMode = failMode(:, iCodeFactor);
        fprintf('Number of IC debonding: %d\n', sum(tmpFailMode==1));
        isIC = tmpFailMode==1;
        %% geometrical properties
        figure;
        subplot(2,2,1);
        tmpEdge = [unique(H_DESIGN_ARRAY_MM)'-1;unique(H_DESIGN_ARRAY_MM)'+1];
        tmpEdge = tmpEdge(:);
        nHBeam = histc(H_DESIGN_ARRAY_MM, tmpEdge);
        nHBeam( nHBeam==0 ) = [];
        nHBeamIc = histc( H_DESIGN_ARRAY_MM(isIC), tmpEdge);
        nHBeamIc( nHBeamIc==0 ) = [];
        bar([nHBeam, nHBeamIc])
        set(gca, 'xticklabel', {num2str(unique(H_DESIGN_ARRAY_MM))})
        xlabel('Beam Height (mm)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')
        
        subplot(2,2,2);
        tmpEdge = [unique(SPAN_DESIGN_ARRAY_MM)'-1;unique(SPAN_DESIGN_ARRAY_MM)'+1];
        tmpEdge = tmpEdge(:);
        nSpan = histc(SPAN_DESIGN_ARRAY_MM, tmpEdge);
        nSpan( nSpan==0 ) = [];
        nSpanIc = histc( SPAN_DESIGN_ARRAY_MM(isIC), tmpEdge);
        nSpanIc( nSpanIc==0 ) = [];
        bar([nSpan, nSpanIc])
        set(gca, 'xticklabel', {num2str(unique(SPAN_DESIGN_ARRAY_MM))})
        xlabel('Beam Span (mm)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')  
        %% concrete properties
        figure;
        tmpEdge = [unique(FC_DESIGN_ARRAY_MPA-FC_DESIGN_ARRAY_MPA*0.01)';
                   unique(FC_DESIGN_ARRAY_MPA+FC_DESIGN_ARRAY_MPA*0.01)'];
        tmpEdge = tmpEdge(:);
        nFc = histc(FC_DESIGN_ARRAY_MPA, tmpEdge); nFc( nFc==0 ) = [];
        nFcIc = histc(FC_DESIGN_ARRAY_MPA(isIC), tmpEdge); nFcIc( nFcIc==0 ) = [];
        bar([nFc, nFcIc])
        set(gca, 'xticklabel', {num2str(unique(FC_DESIGN_ARRAY_MPA))})
        xlabel('Concrete strength (MPa)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')

        %% Steel properties
        figure;
        tmpEdge = [unique(AREA_STEEL_DESIGN_ARRAY_MM2-AREA_STEEL_DESIGN_ARRAY_MM2*0.01)';
                   unique(AREA_STEEL_DESIGN_ARRAY_MM2+AREA_STEEL_DESIGN_ARRAY_MM2*0.01)'];
        tmpEdge = tmpEdge(:);
        nAs = histc(AREA_STEEL_DESIGN_ARRAY_MM2, tmpEdge); nAs( nAs==0 ) = [];
        nAsIc = histc(AREA_STEEL_DESIGN_ARRAY_MM2(isIC), tmpEdge); nAsIc( nAsIc==0 ) = [];
        bar([nAs, nAsIc])
        set(gca, 'xticklabel', {num2str(unique(AREA_STEEL_DESIGN_ARRAY_MM2))})
        xlabel('Steel area (mm^2)')
        ylabel('Number of design cases')
        %% FRP properties
        figure;
        subplot(2,2,1);
        tmpEdge = [unique(E_FRP_DESIGN_ARRAY_MPA-E_FRP_DESIGN_ARRAY_MPA*0.01)';
                   unique(E_FRP_DESIGN_ARRAY_MPA+E_FRP_DESIGN_ARRAY_MPA*0.01)'];
        tmpEdge = tmpEdge(:);
        nEFrp = histc(E_FRP_DESIGN_ARRAY_MPA, tmpEdge); nEFrp( nEFrp==0 ) = [];
        nEFrpIc = histc(E_FRP_DESIGN_ARRAY_MPA(isIC), tmpEdge); nEFrpIc( nEFrpIc==0 ) = [];
        bar([nEFrp, nEFrpIc])
        set(gca, 'xticklabel', {num2str(unique(E_FRP_DESIGN_ARRAY_MPA))})
        xlabel('Modulus of FRP(MPa)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')
        
        subplot(2,2,2);
        tmpEdge = [unique(F_FRP_DESIGN_ARRAY_MPA-F_FRP_DESIGN_ARRAY_MPA*0.01)';
                   unique(F_FRP_DESIGN_ARRAY_MPA+F_FRP_DESIGN_ARRAY_MPA*0.01)'];
        tmpEdge = tmpEdge(:);
        nFFrp = histc(F_FRP_DESIGN_ARRAY_MPA, tmpEdge); nFFrp( nFFrp==0 ) = [];
        nFFrpIc = histc(F_FRP_DESIGN_ARRAY_MPA(isIC), tmpEdge); nFFrpIc( nFFrpIc==0 ) = [];
        bar([nFFrp, nFFrpIc])
        set(gca, 'xticklabel', {num2str(unique(F_FRP_DESIGN_ARRAY_MPA))})
        xlabel('Strength of FRP(MPa)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')
        
        subplot(2,2,3);
        tmpEdge = [unique(T_FRP_DESIGN_ARRAY_MM-T_FRP_DESIGN_ARRAY_MM*0.01)';
                   unique(T_FRP_DESIGN_ARRAY_MM+T_FRP_DESIGN_ARRAY_MM*0.01)'];
        tmpEdge = tmpEdge(:);
        nTFrp = histc(T_FRP_DESIGN_ARRAY_MM, tmpEdge); nTFrp( nTFrp==0 ) = [];
        nTFrpIc = histc(T_FRP_DESIGN_ARRAY_MM(isIC), tmpEdge); nTFrpIc( nTFrpIc==0 ) = [];
        bar([nTFrp, nTFrpIc])
        set(gca, 'xticklabel', {num2str(unique(T_FRP_DESIGN_ARRAY_MM))})
        xlabel('Thickness of FRP (mm)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')
        
        subplot(2,2,4);
        tmpEdge = [unique(B_FRP_DESIGN_ARRAY_MM-B_FRP_DESIGN_ARRAY_MM*0.01)';
                   unique(B_FRP_DESIGN_ARRAY_MM+B_FRP_DESIGN_ARRAY_MM*0.01)'];
        tmpEdge = tmpEdge(:);
        nBFrp = histc(B_FRP_DESIGN_ARRAY_MM, tmpEdge); nBFrp( nBFrp==0 ) = [];
        nBFrpIc = histc(B_FRP_DESIGN_ARRAY_MM(isIC), tmpEdge); nBFrpIc( nBFrpIc==0 ) = [];
        bar([nBFrp, nBFrpIc])
        set(gca, 'xticklabel', {num2str(unique(B_FRP_DESIGN_ARRAY_MM))})
        xlabel('Width of FRP (mm)')
        ylabel('Number of design cases')
%         legend('Original design cases', 'IC debonding')
    otherwise
end