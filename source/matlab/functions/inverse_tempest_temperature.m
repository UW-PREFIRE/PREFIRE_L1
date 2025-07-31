%% TEMPEST-D Temperature Conversion

function [itemp] = inverse_tempest_temperature(tdouble)

B=3891; %coefficients
R0=5000;    %resistor 0
T0=298.15;  %zeroth temperature

rinf=R0*exp(-B/T0);
res=rinf*(exp(B./tdouble));
iread=5./(res+10000);
vread=-(iread.*10000.-5);
itemp=round(4096*vread./5);

end

