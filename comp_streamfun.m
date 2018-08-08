function [PSI,psi_min,psi_max] = comp_streamfun(q,Grid)

%% Q1: Part 3

% author: Evan Ramos
% date:   26 October 2015

% Description: The streamfunction of a numerical solution is computed at
% the cell corners given the vector of numerical fluxes across
% the cell faces.

% Input: q = [qx;qy] Nf by 1 vector of all face fluxes. x-fluxes first.
% Grid = structure containing all information about the grid.

% Output: PSI = Ny+1 by Nx+1 matrix containing the values of the
% streamfunction along the cell corners.
% psi_min = minimum of PSI
% psi_max = maximum of PSI

if strcmp(Grid.psi_x0,'xmin_ymin') && strcmp(Grid.psi_dir,'xy')
    % Integrate horizontally along ymin boundary first then vertically into the domain.
    Qymin = [0 cumsum(Grid.A(Grid.dof_f_ymin).*q(Grid.dof_f_ymin))']; % Integral of flow into 
    Qx = reshape(Grid.A(1:Grid.Nfx).*q(1:Grid.Nfx),Grid.Ny,Grid.Nx+1); % Horizontal fluxes
    PSI = cumsum([-Qymin;Qx],1); % integrals into domain with Qymin as initial value
elseif strcmp(Grid.psi_x0,'xmin_ymin') && strcmp(Grid.psi_dir,'yx')
    % Integrate vertically along xmin boundary first then horizontally into the domain.
    Qxmin = [0; cumsum(Grid.A(Grid.dof_f_xmin).*q(Grid.dof_f_xmin))]; % Integral of flow into xmin boundary
    Qy = reshape(Grid.A(Grid.Nfx + (1:Grid.Nfy)).*q(Grid.Nfx+(1:Grid.Nfy)),Grid.Ny+1,Grid.Nx); % Vertical fluxes
    PSI = cumsum([-Qxmin Qy],2); % integrals into domain with Qxmin as initial value
else
    error('Unknown streamfunction setup.\n')
end

psi_min = min(PSI(:)); psi_max = max(PSI(:));

end