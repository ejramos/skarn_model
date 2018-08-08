function [Rho,  Mu] = steam_table(T,P,Grid)
%%
%author: Evan J. Ramos
%date:   22 Jan 2017

%Description: This function loads steam table data regarding fluid density,
%and dynamic viscosity and given a set of temperatures and
%pressures, returns interpolated values from matrices.

%Inputs:
%         T    --> Grid.N x 1 vector containing temperatures [K]
%         P    --> Grid.N x 1 vector containing pressures    [Pa]
%         Grid --> structure containing the grid properties

%Outputs:
%         Rho  --> Grid.N x 1 vector containing fluid density [kg/m^3]
%         Cp   --> Grid.N x 1 vector containing heat capacity [J/mol/K]
%         Mu   --> Grid.N x 1 vector containing viscosity     [Pa*s]
%

%% Load steam tables

% NOTE:
% Density: from Burnham et al. (1969)
% Viscosity: from Haar et al. (1984)

density   = load('steam_table_density.txt');
viscosity = load('steam_table_viscosity.txt')';

%% Make interpolated matrices from steam table data

%density
t_rho                    = density(2:end,1);
p_rho                    = density(1,2:end);
[P_rho, T_rho]           = meshgrid(p_rho,t_rho);

t_rho_fine               = linspace(min(t_rho),max(t_rho));
p_rho_fine               = linspace(min(p_rho),max(p_rho));
[P_rho_fine, T_rho_fine] = meshgrid(p_rho_fine,t_rho_fine);

rho_fine                 = interp2(P_rho,T_rho,density(2:end,2:end),...
                                   P_rho_fine,T_rho_fine,'spline');
rho_fine(rho_fine <= 400) = 400;

%viscosity
t_vis                    = viscosity(2:end,1);
p_vis                    = viscosity(1,2:end);
[P_vis, T_vis]           = meshgrid(p_vis,t_vis);

t_vis_fine               = linspace(min(t_vis),max(t_vis));
p_vis_fine               = linspace(min(p_vis),max(p_vis));
[P_vis_fine, T_vis_fine] = meshgrid(p_vis_fine,t_vis_fine);

vis_fine                 = interp2(P_vis,T_vis,viscosity(2:end,2:end),...
                                   P_vis_fine,T_vis_fine,'spline');

%% Locate temperature and pressure in fine vectors --> new fluid properties

Rho = zeros(Grid.N,1);
Mu  = zeros(Grid.N,1);

for i = 1:Grid.N
    %density
    [~, r_rho] = min(abs(T(i)-t_rho_fine));
    [~, c_rho] = min(abs(P(i)-p_rho_fine));
    Rho(i)     = rho_fine(r_rho,c_rho);
    
    %viscosity
    [~, r_vis] = min(abs(T(i)-t_vis_fine));
    [~, c_vis] = min(abs(P(i)-p_vis_fine));
    Mu(i)      = vis_fine(r_vis,c_vis);
end

end