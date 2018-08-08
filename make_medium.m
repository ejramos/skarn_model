function [Mineral, phi, Kd, varargout] = make_medium(filename,Grid,sett,pluton)
%% File description

%author: Evan J. Ramos
%date:   12 Jan 2017

%This function creates all variables associated with medium properties,
%including porosity, permeability (which is a function of porosity in all
%areas except those with faults), and the distribution of minerals and
%their relative abundances.

%Inputs:
%filename --> string of the name of the file contains all known properties
%             of the minerals used in the simulation
%Grid     --> structure containing all properties of the grid
%set      --> an integer which represents an option for the "setting" in
%             which would like to be modeled.
%pluton   --> an integer which represents an option for the placement of a
%             pluton spatially, either on the bottom left, bottom middle,
%             or bottom right of the domain. Pluton size solely depends
%             upon the grid size


%Outputs:
%Mineral  --> structure containing all mineral properties found within the
%             file passed to the function
%phi      --> Grid.N x 1 vector containing all of the porosity values for
%             the grid cells
%Kd       --> Grid.Nf x Grid.Nf sparse diagonal matrix containing the
%             harmonic mean of the permeability values at cell faces
%Pluton   --> structure containing all of the degrees of freedom associated
%             with the pluton

%Example call:
% >> [Mineral, phi, Kd] = make_medium('mineralisotopecontsant.txt',...
%                                     Grid,3,0);
%     "Given the mineral properties found in mineralisotopeconstant.txt
%     and the grid properties found in the structure Grid, generate a
%     medium that follows setting option 3 and
%     pluton location option 0 (no pluton)."

% The way that the porosities and permeabilities are implemented follow the
% way in which these structures develop. First there is the bulk crustal
% porosity and permeability, representative of the pre-existing basin in
% which fluid is flowing prior to skarn formation. Then, the pluton is 
% emplaced and overprints any of the structures that exist.

%% Build mineral structure

Mineral             = build_mineral_new(filename,Grid);

%% Determine porosity field

phi_min = 1e-4;                       %[m^3/m^3] minimum imposed porosity
phi_0   = 5e-3;                       %[m^3/m^3] reference porosity
k_0     = 1e-17;                      %[m^2] reference permeability

n       = 3;                          %integer exponent for power law calc

phi_p   = 2.5e-3;                     %[m^3/m^3] pluton porosity
phi_lim = 1e-2;                       %[m^3/m^3] limestone porosity

switch sett
    case 1
        %Uniform limestone country rock --> uniform permeability
        phi        = phi_lim*ones(Grid.N,1);
        k          = perm_from_poro(k_0,phi,phi_min,phi_0,n);
    case 2
        % permeability decrease of 3 orders of magnitude with depth
        
        depth             = flipud(Grid.yc); 
        k                 = logspace(-16,-13,length(depth));
        k                 = reshape(repmat(k,1,Grid.Nx),Grid.N,1);
        phi               = poro_from_perm(k_0,k,phi_min,phi_0,n);
        
end

%% Determine pluton location

switch pluton
    case 1
        [dof_pluton,dof_f_pluton,...
            bound, dof_rock] = make_pluton(Grid,'left');
    case 2
        [dof_pluton,dof_f_pluton,...
            bound,dof_rock] = make_pluton(Grid,'center');
    case 3
        [dof_pluton,dof_f_pluton,...
            bound,dof_rock] = make_pluton(Grid,'right');
end


if pluton ~= 0
    varargout{1}.dof       = dof_pluton;
    varargout{1}.dof_f     = dof_f_pluton;
    varargout{1}.dof_rock  = dof_rock;
    varargout{1}.border    = bound;
    phi(varargout{1}.dof)  = phi_p;
    k(varargout{1}.dof)    = perm_from_poro(k_0,phi_p,phi_min,phi_0,n);
end

%% Generate permeability matrix

Kd = comp_mean(reshape(k,Grid.Ny,Grid.Nx),-1,1,Grid);
if nargout == 5
    varargout{2} = reshape(k,Grid.Ny,Grid.Nx);
end

%% Create Phi field for Mineral structure

for j = 1:Grid.N
    %pluton mineral compositions
    if exist('dof_pluton','var')
        if ismember(j,dof_pluton)
            Mineral.cal.phi(j) = 0;
            Mineral.dol.phi(j) = 0;
            Mineral.qtz.phi(j) = .35;
            Mineral.An40.phi(j) = .65;
            Mineral.mus.phi(j) = 0;
            Mineral.bte.phi(j) = 0;
            Mineral.hbd.phi(j) = 0;
        end
    end
    %if it's not the location of a pluton, it'll be the
    %country rock (limestone)
    if isequal(Mineral.cal.phi(j),-1)
        Mineral.cal.phi(j) = .75;
        Mineral.dol.phi(j) = .25;
        Mineral.qtz.phi(j) = 0;
        Mineral.An40.phi(j) = 0;
        Mineral.mus.phi(j) = 0;
        Mineral.bte.phi(j) = 0;
        Mineral.hbd.phi(j) = 0;
    end

end
end