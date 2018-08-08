function phi = poro_from_perm(k_0,k,phi_min,phi_0,n)

%author: Evan J. Ramos
%date:   12 Jan 2017

%Description: This function calculates crustal porosity based upon a
%             power law for permeabilty (following Walder and Nur, 1984).

%Inputs:
%
% k_0      --> reference permeability value [m^2]
% k        --> input permeability value [m^2]
% phi_min  --> minimum porosity value [m^3/m^3]
% phi_0    --> reference porosity value [m^3/m^3]
% n        --> integer value for the power law
%              a value of 2: representative of channel flow
%              a value of 3: representative of parallel fracture evolution

%Outputs:
% phi      --> computed porosity [m^3/m^3]

phi   = ((k/k_0).*(phi_0.^n - phi_min.^n) + phi_min.^n).^(1/n);

end