function [q] = comp_flux_p(D,Kd,zhat,G,p,fs,Grid,Param)
% author: Evan J. Ramos, Jake S. Jordan
% date:   08 June 2016
% Description:
% Computes the mass conservative fluxes across all boundaries from the 
% residual of the compatability condition over the boundary cells.
% Note: Current implmentation works for all cases where one face 
%       is assigend to each bnd cell. So conrner cells must have
%       natural BC's on all but one face.
%
% Input:
% D = N by Nf discrete divergence matrix.
% Kd = Nf by Nf conductivity matrix.
% G = Nf by N discrete gradient matrix.
% vert = N by 1 vector of depths.
% p = N by 1 vector of fluid pressure in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% Param = structure contaning problem paramters and information about BC's
%
% Output:
% q = Nf by 1 vector of fluxes across all cell faces

%% Compute interior fluxes

if isempty(zhat)
    zhat = zeros(Grid.Nf,1);
end
q = -Kd*(G*p + zhat);

%rewrite code such that you pass the function q to the comp_flux equation
%in the form of a function HANDLE. 

%% Compute boundary fluxes
dof_cell = [Param.dof_dir;Param.dof_neu];
dof_face = [Param.dof_f_dir;Param.dof_f_neu];
sign = ismember(dof_face,[Grid.dof_f_xmin;Grid.dof_f_ymin])...
       -ismember(dof_face,[Grid.dof_f_xmax;Grid.dof_f_ymax]);

if ~isempty(dof_cell)
    q(dof_face) =  sign.*( D(dof_cell,:) * q - fs(dof_cell)).*Grid.V(dof_cell)./Grid.A(dof_face);
end
end