function [D,G,I,varargout]=build_ops(Grid) 

%% Q3

% author:  Evan J. Ramos
% date:    19 October 2015 (most recently updated for 2D compatibility)
% updated: 02 June 2016

% Description: 
% This function computes the discrete divergence and gradient matrices on a 
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions. 

% Input: 
% Grid = structure containing all pertinent information about the grid. 

% Output: 
% D = Nx by Nx+1 discrete divergence matrix (if 1D)

B = [-1*ones(Grid.Nx,1) ones(Grid.Nx,1)]; %diagonal elements
D = (1/Grid.dx)*spdiags(B,[0 1],Grid.Nx,Grid.Nx+1);

% in the case that the flow is two-dimensional
if ~isfield(Grid,'D2')
    
    if Grid.Nx > Grid.Ny
        D_part = (Grid.dx/Grid.dy)*D(1:Grid.Ny,1:(Grid.Ny+1));
    else
        B = [-1*ones(Grid.Ny,1) ones(Grid.Ny,1)];
        D_part = (1/Grid.dy)*spdiags(B,[0 1],Grid.Ny,Grid.Ny+1);
    end
    
    Dy = spblkdiag(D_part,Grid.Nx);
    B = [-1*ones(Grid.N,1) ones(Grid.N,1)]; %diagonal elements
    Dx = (1/Grid.dx)*spdiags(B,[0 Grid.Ny],Grid.N,Grid.Nfx);
    D = [Dx Dy];

end
    
% G = Nx+1 by Nx discrete gradient matrix
% Gx = Nf by N matrix containing just Gx matrix --> rest are zeros (if 2D)
% Gy = Nf by N matrix containing just Gy matrix --> rest are zeros (if 2D)
G = -D';

if isfield(Grid,'D2')
    G(1) = 0;
    G(end) = 0;
else
    bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax; ...
        Grid.dof_f_ymin; Grid.dof_f_ymax]';
    G(bnd,:) = 0;
end

    if nargout == 5
        %Gx
        varargout{1} = G;
        varargout{1}(Grid.Nfx+1:end,:) = 0;
        %Gy
        varargout{2} = G;
        varargout{2}(1:Grid.Nfx,:) = 0;
    end

% I = N by N identity matrix 
I = speye(Grid.N);

end