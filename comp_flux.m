function [q] = comp_flux(D,Kd,G,h,fs,Grid,Param)
% author: Marc Hesse
% date: 25 Nov 2014, 10 Jul 2015
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
% h = N by 1 vector of flow potential in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% Param = structure contaning problem paramters and information about BC's
%
% Output:
% q = Nf by 1 vector of fluxes across all cell faces
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> Param.dof_dir = Grid.dof_xmin;
% >> Param.dof_f_dir = Grid.dof_f_xmin;
% >> g = 0;
% >> Param.dof_neu = []; Param.dof_f_neu =[];
% >> [B,N,fn] = build_bnd(Param,Grid);
% >> h = solve_lbvp(L,fs+fn,B,g,N);
% >> q = comp_flux(D,1,G,h,fs,Grid,Param);

%% Compute interior fluxes
q = -Kd*G*h;

if isfield(Param,'flow')
    Param.dof_dir = Param.flow.dof_dir;
    Param.dof_f_dir = Param.flow.dof_f_dir;
    Param.dof_neu = Param.flow.dof_neu;
    Param.dof_f_neu = Param.flow.dof_f_neu;
elseif isfield(Param,'trans')
    Param.dof_dir = Param.trans.dof_dir;
    Param.dof_f_dir = Param.trans.dof_f_dir;
    Param.dof_neu = Param.trans.dof_neu;
    Param.dof_f_neu = Param.trans.dof_f_neu;
end

%% Compute boundary fluxes
dof_cell = [Param.dof_dir;Param.dof_neu];
dof_face = [Param.dof_f_dir;Param.dof_f_neu];
sign = ismember(dof_face,[Grid.dof_f_xmin;Grid.dof_f_ymin])...
      -ismember(dof_face,[Grid.dof_f_xmax;Grid.dof_f_ymax]);

q(dof_face) =  sign.*( D(dof_cell,:) * q - fs(dof_cell)).*Grid.V(dof_cell)./Grid.A(dof_face);
end