function [A] = flux_upwind(q,Grid)
% author: Marc Hesse
% date: 15 April 2015
% Description:
% This function computes the upwind flux matrix from the flux vector.
%
% Input:
% q = Nf by 1 flux vector from the flow problem.
% Grid = structure containing all pertinent information about the grid.
%
% Output:
% A = Nf by Nf matrix contining the upwinded fluxes
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> q = ones(Grid.Nf,1);
% >> [A] = flux_upwind(q,Grid);

Nx = Grid.Nx; Ny = Grid.Ny; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces

if ((Nx>1) && (Ny==1)) || ((Nx==1) && (Ny>1)) % 1D
    %% One dimensinal
    % too make this work for 1D in y-dir need to replace Nx with N!
    qn = min(q(1:Nx),0);
    qp = max(q(2:Nx+1),0);
    A = spdiags([qp,qn],[-1 0],Grid.Nx+1,Grid.Nx);
elseif (Nx>1) && (Ny>1) % 2D
    %% Two dimensional
    % x - fluxes
    qxn = min(q(1:Nfx-Ny),0); 
    qxp = max(q(Ny+1:Nfx),0);
    Ax = spdiags([qxp,qxn],[-Ny 0],Nfx,N);
    
    % y-fluxes
    QY = reshape(q(Nfx+1:end),Grid.Ny+1,Grid.Nx);
    qyn = min(reshape(QY(1:Grid.Ny,:),Grid.N,1),0);
    qyp = max(reshape(QY(2:Grid.Ny+1,:),Grid.N,1),0);
    row_p = (Grid.Ny+1)*repmat([0:Grid.Nx-1],Grid.Ny,1)+repmat([2:Grid.Ny+1]',1,Grid.Nx); 
    row_n = (Grid.Ny+1)*repmat([0:Grid.Nx-1],Grid.Ny,1)+repmat([1:Grid.Ny]',  1,Grid.Nx); 
    Ay = sparse([row_p(:);row_n(:)],[Grid.dof';Grid.dof'],[qyp;qyn]);
    
    A = [Ax; Ay];
end