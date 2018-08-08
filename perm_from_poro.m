function k = perm_from_poro(k_0,phi,phi_min,phi_0,n)

%author: Evan J. Ramos
%date:   12 Jan 2017

%Description: This function calculates crustal permeability based upon a
%             power law for porosity (following Walder and Nur, 1984).

%Inputs:
%
% k_0      --> reference permeability value [m^2]
% phi      --> input porosity value [m^3/m^3]
% phi_min  --> minimum porosity value [m^3/m^3]
% phi_0    --> reference porosity value [m^3/m^3]
% n        --> integer value for the power law
%              a value of 2: representative of channel flow
%              a value of 3: representative of parallel fracture evolution

%Outputs:
% k        --> computed permeability [m^2]

k                 = k_0.*((phi.^n - phi_min.^n)./(phi_0.^n - phi_min.^n));

k(phi <= phi_min) = 0;

end