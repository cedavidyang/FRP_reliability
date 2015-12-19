function beta = life2re(life, beta_T_50)
% return design service life based on reliabiltiy index beta

beta = norminv((normcdf(beta_T_50)).^(50./life));

return
end