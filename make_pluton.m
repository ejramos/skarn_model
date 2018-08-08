function [dir_locs,neu_locs,boundary,dof_rock] = make_pluton(Grid,type)

% author: Evan J. Ramos
% date:   08 May 2016

% Description: This function creates a permeability matrix of ones and
% zeros where ones represent the location of a pluton and zeros are the
% rest of the domain. It also returns the Neumann boundary locations (for
% constant heat fluxes) at the pluton country rock interface.
%
% The pluton is rectangular in shape

% Inputs:   Grid --> structure containing properties of the grid
%           type --> string denoting location of pluton

% Outputs:
%           dir_locs     --> length(bnd) x 1 vector containing the
%                            degrees of freedom of the cell centers
%                            associated with the entire pluton
%           neumann_locs --> length(dof_f_neu) x 1 vector containing the
%                            degrees of freedom of the cell faces 
%                            associated with the outermost pluton body
%           boundary     --> length(bnd) x 1 vector containing the
%                            degrees of freedom of the cell centers
%                            associated with the outermost pluton body
%           dof_rock     --> matrix of (x,y) coords directly above and
%                            beside the pluton (used for tracking how delta
%                            values vary through time at different points
%                            in space)

cell_center_mat     = reshape(Grid.dof,Grid.Ny,Grid.Nx);

%% Defining pluton dimensionality based on cell location

row_pluton          = floor(.25*Grid.xmax);  %width of pluton (x-distance)
col_pluton          = floor(.40*Grid.ymax);  %length of pluton(y-distance)

%{xmin,xmax,ymin,ymax} --> row, column location of pluton
ymin                = 1;
[~,ymax]            = min(abs(Grid.yc-col_pluton));

if strcmpi(type,'left')
    
    xmin            = 1;
    [~,xmax]        = min(abs(Grid.xc-row_pluton));
    rocks           = zeros(4,2);
    rocks(1,:)      = [ymax + floor(.04*Grid.dy), xmin + ceil(.3*(xmax-xmin))];
    rocks(2,:)      = [ymax + floor(.04*Grid.dy), xmax - ceil(.3*(xmax-xmin))];
    rocks(3,:)      = [floor(.3*(ymax-ymin)), xmax + floor(.04*Grid.Nx)];
    rocks(4,:)      = [floor(.7*(ymax-ymin)), xmax + floor(.04*Grid.Nx)];
 
elseif strcmpi(type,'right')
    
    row_pluton      = Grid.xmax - row_pluton;
   
    [~,xmin]        = min(abs(Grid.xc-row_pluton)); 
    xmax            = Grid.Nx;
    rocks           = zeros(4,2);
    rocks(1,:)      = [ymax + floor(.04*Grid.dy), xmin + ceil(.3*(xmax-xmin))];
    rocks(2,:)      = [ymax + floor(.04*Grid.dy), xmax - ceil(.3*(xmax-xmin))];
    rocks(3,:)      = [floor(.3*(ymax-ymin)), xmin - ceil(.04*Grid.Nx)];
    rocks(4,:)      = [floor(.7*(ymax-ymin)), xmin - ceil(.04*Grid.Nx)];
 
    
elseif strcmpi(type,'center')
    
    row_pluton_min  = floor(Grid.xmax/2 - row_pluton/2);
    row_pluton_max  = floor(Grid.xmax/2 + row_pluton/2);
    
    [~,xmin]        = min(abs(Grid.xc-row_pluton_min)); 
    [~,xmax]        = min(abs(Grid.xc-row_pluton_max));
    rocks           = zeros(6,2);
    rocks(1,:)      = [ymax+1, xmin + ceil(.3*(xmax-xmin))];
    rocks(2,:)      = [ymax+1, xmax - ceil(.3*(xmax-xmin))];
    rocks(3,:)      = [floor(.3*(ymax-ymin)), xmax+2];
    rocks(4,:)      = [floor(.7*(ymax-ymin)), xmax+2];
    rocks(5,:)      = [floor(.3*(ymax-ymin)), xmin-2];
    rocks(6,:)      = [floor(.7*(ymax-ymin)), xmin-2];
end
     
%% Rock locations
dof_rock = zeros(length(rocks),1);
for i = 1:length(dof_rock)
    dof_rock(i) = cell_center_mat(rocks(i,1),rocks(i,2));
end
%% Neumann locations

% %Get cell locations of the boundaries
top                 = cell_center_mat(ymax,xmin:xmax);
bottom              = cell_center_mat(ymin,xmin:xmax);
left                = cell_center_mat(ymin:ymax,xmin);
right               = cell_center_mat(ymin:ymax,xmax);
% 
% %unique --> removes repeats
% %excluding cells that are on the border of the 2D domain
% 
boundary             = unique([top'; bottom'; left; right]);
boundary             = boundary(~ismember(boundary,Grid.dof_xmin));
boundary             = boundary(~ismember(boundary,Grid.dof_xmax));
boundary             = boundary(~ismember(boundary,Grid.dof_ymin));

dir_locs             = cell_center_mat(ymin:ymax,xmin:xmax);
dir_locs             = reshape(dir_locs,numel(dir_locs),1);

%get degrees of freedom --> reshape
dof_fs               = get_dof_faces2D(dir_locs,Grid);
neu_locs             = reshape(dof_fs,numel(dof_fs),1);

end