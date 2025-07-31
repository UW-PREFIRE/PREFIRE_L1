function [ output] = exitance( lam1, lam2, temp )
%This function integrates the spectral radiance of the Planck function
%using the approximation in http://www.spectralcalc.com/blackbody/inband_radiance.html
%input units are um, um, K
a_1 = 8.73068E-13;
%output consistent with W/cm^2, but it is not 'per micron' scaled, 
%unless lam2-lam1 = 1.0 um
c_2 = 14387.69;  % (1/k) Boltzmann constant in um K  10000/.695 cm-1 K-1
    sum1 = 0;
    sum2 = 0;
    
   
       x = c_2 / lam1 / temp;
       for N=1:100
           y = N * x;
           term = (exp(-y) / N ^ 4) * (6 + y * (6 + y * (3 + y)));
           sum1 = sum1 + term;
           % 'If term / sum1 < TOL Then Exit For
       end
   
   

       x = c_2 / lam2 / temp;
       for N = 1:100
           y = N * x;
           term = (exp(-y) / N ^ 4) * (6 + y * (6 + y * (3 + y)));
           sum2 = sum2 + term;
           %'If term / sum2 < TOL Then Exit For
       end


    
    output = a_1 * temp ^ 4 * (sum2 - sum1);
end

