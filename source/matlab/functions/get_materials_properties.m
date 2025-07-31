function [materials] = get_materials_properties()

materials.R_o=0.95;            % this is 1 - emissivity of polished aluminum
materials.R_f=0.95;            % estimate for filter reflectivity

%these are materials properties
materials.em = 0.77; %emissivity of blocking filter (same as anodized Al for now)
materials.R_C = 0.3; % diamond reflectance
materials.Z_e = 0.2; % zero order grating efficiency
materials.e_Al = 0.77; % anodized (rough is 0.07, oxidized is 0.25, polished is 0.04
materials.e_C = 0.15; % good assumption for longer wavelengths, smaller bigger for thicker material
materials.e_MB = 0.99; % emissivity of Martin Black paint (needed to bring down emittance from outFOV)
materials.e_GB = 0.99; % emissivity of gold black on detectors

end
