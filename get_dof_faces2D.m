function dof_fs = get_dof_faces2D(dof,Grid)

% author: Evan J. Ramos
% date:   07 May 2016

% Description: This function takes in the degrees of freedom associated
% with the cell center and the Grid parameters in order to return the
% degrees of freedom associated with the cell faces.

% Inputs:  dof    --> number of the cell center
%          Grid   --> structure containing properties of the grid

% Outputs: dof_fs --> 4 x length(dof) vector containing the degrees of 
%                     freedom associated with the cell faces
%                     
%                     Elements 1 & 2: X faces
%                     Elements 3 & 4: Y faces

dof_fs      = zeros(4,length(dof));

dof_fs(1,:) = dof;
dof_fs(2,:) = dof + Grid.Ny;

cell_mat    = flipud(reshape(Grid.dof,Grid.Ny,Grid.Nx));
for i = 1:length(dof)
    
    [~,c]       = find(dof(i) == cell_mat);
    dof_fs(3,i) = Grid.Nfx + dof(i) + (c-1);
    dof_fs(4,i) = dof_fs(3,i) + 1;

end

end

