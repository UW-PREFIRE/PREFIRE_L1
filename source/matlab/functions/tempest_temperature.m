%% TEMPEST-D Temperature Conversion

function [ttemp] = tempest_temperature(tdata)

B=3891; %coefficients
R0=5000;    %resistor 0
T0=298.15;  %zeroth temperature

vread=5.*tdata./4096;
iread=(5-vread)./10000;
res=(5./iread)-10000;
rinf=R0*exp(-B/T0);
ttemp=B./log(res./rinf);

end
