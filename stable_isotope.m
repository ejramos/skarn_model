function [alpha, delta_out] = stable_isotope(D,E,F,T,delta_in,numden)

% author: Evan J. Ramos
% date:   12 January 2016

% Description:
% This function computes the stable isotope value (delta 18-O in
% particular) based upon a temperature T and the constants associated with
% the particular phase relationships (D,E,F). The fractionation factor
% (alpha) is calculated from a ratio of one delta value with another from
% the differing reservoirs. The input 'numden' is a string which indicates 
% whether the input delta value pertains to the numerator or the 
% denominator in the fractionation factor for the given calculation.

% T [K]
% delta [per-mil]
% D,E,F [unitless]


%%%%%%% Assumptions %%%%%%%
% 1. Although involving external fluids, I'm assuming closed system
% processes. Thus, equilibrium is assumed.


alpha = exp(((D*10^6./(T.^2)) + (E*10^3./T) + F)/1000);


if strcmp(numden,'num')
    delta_out = (1000+delta_in)./alpha - 1000;
elseif strcmp(numden,'den')
    delta_out = (1000+delta_in).*alpha - 1000;
end

end