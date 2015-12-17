function [re_smp, mean_RE, std_RE, upper_RE, lower_RE, norm_RE_yr30, ...
    norm_RE_yr40, norm_RE_yr50] = get_re_data( data_struct, psi_f, code, beta_T_50)
% get re and summaries from data
npsi = length(psi_f);
beta_T_30 = life2re(30, beta_T_50);
beta_T_40 = life2re(40, beta_T_50);
switch lower(code)
    case {'aci', 'acinew'}
        for ipsi = 1:npsi
            re_smp = data_struct.re_data{1}(:, ipsi, :);
            re_smp = re_smp(:);
            re_smp( isnan(re_smp) ) = [];
            norm_RE_yr30(ipsi) = mean((re_smp-beta_T_30).^2);
            norm_RE_yr40(ipsi) = mean((re_smp-beta_T_40).^2);
            norm_RE_yr50(ipsi) = mean((re_smp-beta_T_50).^2);
            mean_RE(ipsi) = mean(re_smp);
            std_RE(ipsi) = std(re_smp);
            upper_RE(ipsi) = max(re_smp);
            lower_RE(ipsi) = min(re_smp);
        end
    otherwise
        for ipsi = 1:npsi
            try
                re_smp = data_struct.re_data{ipsi}(:, 1, :);
            catch
                re_smp = data_struct.re_data(ipsi, :, :, :);
            end
            re_smp = re_smp(:);
            re_smp( isnan(re_smp) ) = [];
            norm_RE_yr30(ipsi) = mean((re_smp-beta_T_30).^2);
            norm_RE_yr40(ipsi) = mean((re_smp-beta_T_40).^2);
            norm_RE_yr50(ipsi) = mean((re_smp-beta_T_50).^2);
            mean_RE(ipsi) = mean(re_smp);
            std_RE(ipsi) = std(re_smp);
            upper_RE(ipsi) = max(re_smp);
            lower_RE(ipsi) = min(re_smp);
        end
end

return
end
