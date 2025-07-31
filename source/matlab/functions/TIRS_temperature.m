function [T_val] = TIRS_temperature(DN)
% Convert TIRS thermistor DN(s) to temperature value(s)

R_pullup = 10.;  % [kohm]
R_0 = 5.;  % [kohm]
T_0_inv = 1./298.15;  % [1/K]
beta_inv = 1./3891.;  % [1/K]

DN_term = log((R_pullup/R_0).*(DN./(4096-DN)));  % [-]
T_inv = T_0_inv+beta_inv.*DN_term;  % [1/K]

T_val = 1./T_inv;  % [K]

end
