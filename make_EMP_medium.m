function [Mineral, phi, Kd, k, Pluton] = make_EMP_medium(filename,Grid,Emp_coords)
%% File description

%author: Evan J. Ramos
%date:   09 Jan 2018

%This function creates all variables associated with medium properties,
%including porosity, permeability (which is a function of porosity in all
%areas except those with faults), and the distribution of minerals and
%their relative abundances for the newest formulation of Empire Mountain.

%Inputs:
%filename   --> string of the name of the file contains all known
%               properties of the minerals used in the simulation
%Grid       --> structure containing all properties of the grid
%Emp_coords --> structure containing the coordinate locations for each unit
%               in the Empire Mountain set up

%Outputs:
%Mineral    --> structure containing all mineral properties found within
%               the file passed to the function
%phi        --> Grid.N x 1 vector containing all of the porosity values for
%               the grid cells
%Kd         --> Grid.Nf x Grid.Nf sparse diagonal matrix containing the
%               harmonic mean of the permeability values at cell faces
%k          --> Grid.N x 1 vector containing the permeability of each cell
%Pluton     --> structure containing all of the degrees of freedom
%               associated with the pluton

%% Build mineral structure
Mineral             = build_mineral_new(filename,Grid); %main variable
Mineral.unit        = cell(Grid.N,1); %adding cell array with unit classif.
[Xc, Yc]            = meshgrid(Grid.xc,Grid.yc);

%% Porosity constraints

%Walder and Nur (1984) relationships
phi_min = 1e-4;                       %[m^3/m^3] minimum imposed porosity
phi_0   = 5e-3;                       %[m^3/m^3] reference porosity
k_0     = 1e-17;                      %[m^2] reference permeability
n       = 3;                          %integer exponent for power law calc
n_fault = 2;                          %integer exponent for fault calc

%Feature porosity
phi_p   = 2.5e-3;                     %[m^3/m^3] pluton porosity
phi_f   = 5e-2;                       %[m^3/m^3] fault porosity
%phi_lim = .01;                       %[m^3/m^3] limestone porosity

%Rock type-dependent porosity relationships
                     %Approximated from Allen and Allen, 2013 p. 333 vals

%surface porosity
phi_surf.cs   = .3;  %[m^3/m^3] calc-silicate surface porosity
phi_surf.m    = .4;  %[m^3/m^3] carbonate surface porosity
phi_surf.ig   = .25; %[m^3/m^3] extrusive (rhy,dac,and) surface porosity
phi_surf.far  = .60; %[m^3/m^3] Farewell unit (phyllite/hornfels) surface porosity

%porosity depth coefficient
phi_dcof.cs   = .4;  %[km^-1] calc-silicate poro depth coeff
phi_dcof.m    = .6;  %[km^-1] carbonate poro depth coeff
phi_dcof.ig   = .45; %[km^-1] extrusive poro depth coeff
phi_dcof.far  = .5;  %[km^-1] Farewell unit poro depth coeff

%% 1. Create uniform limestone medium over entire domain

%First establish the modal mineralogy: pure limestone
Mineral.cal.phi = .75*ones(Grid.N,1);
Mineral.dol.phi = .25*ones(Grid.N,1);
Mineral.qtz.phi = zeros(Grid.N,1);
Mineral.An40.phi = zeros(Grid.N,1);
Mineral.mus.phi = zeros(Grid.N,1);
Mineral.bte.phi = zeros(Grid.N,1);
Mineral.hbd.phi = zeros(Grid.N,1);

%Second establish porosity (based off Allen and Allen, 2013)
k      = reshape(repmat(logspace(-16,-13,Grid.Ny),1,Grid.Nx),Grid.N,1);
phi    = poro_from_perm(k_0,k,phi_min,phi_0,n);

%Third classify each cell as limestone
Mineral.unit(true(Grid.N,1)) = {'m'};



%% 2. Edit region i

%Build Farewell unit
%Locate it on the prescribed grid (get its dofs)
dof_temp      = Grid.dof(...
                (Xc > Emp_coords.I(1) & Xc < Emp_coords.I(2)) & ...
                (Yc > Emp_coords.I(3) & Yc < Emp_coords.I(4)));
dof_tot       = length(dof_temp);
d_col         = Grid.yc(Grid.yc > Emp_coords.I(3) & ...
                        Grid.yc < Emp_coords.I(4));
%calculate depth-dependent poro for a column of Farewell
d_col         = (Grid.ymax - d_col)/1000;
phi_d         = phi_surf.far*exp(-phi_dcof.far*d_col);
%replicate it for entire region
phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
phi(dof_temp) = phi_total;
%Update unit classification
Mineral.unit(dof_temp) = {'Farewell'};


%Build Kmrt
%Locate it on the prescribed grid (get its dofs)
dof_temp      = Grid.dof(...
                (Xc > Emp_coords.i.Kmrt(1) & Xc < Emp_coords.i.Kmrt(2)) & ...
                (Yc > Emp_coords.i.Kmrt(3) & Yc < Emp_coords.i.Kmrt(4)));
dof_tot       = length(dof_temp);
d_col         = Grid.yc(Grid.yc > Emp_coords.i.Kmrt(3) & ...
                        Grid.yc < Emp_coords.i.Kmrt(4));
%calculate depth-dependent poro for a column of Kmrt
d_col         = (Grid.ymax - d_col)/1000;
phi_d         = phi_surf.ig*exp(-phi_dcof.ig*d_col);
%replicate it for entire region
phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
phi(dof_temp) = phi_total;
%Update unit classification
Mineral.unit(dof_temp) = {'Kmrt'};



%% 3. Edit region ii

%Build Farewell unit
%Locate it on the prescribed grid (get its dofs)
dof_temp      = Grid.dof(...
                (Xc > Emp_coords.II(1) & Xc < Emp_coords.II(2)) & ...
                (Yc > Emp_coords.II(3) & Yc < Emp_coords.II(4)));
dof_tot       = length(dof_temp);
d_col         = Grid.yc(Grid.yc > Emp_coords.II(3) & ...
                        Grid.yc < Emp_coords.II(4));
%calculate depth-dependent poro for a column of Farewell
d_col         = (Grid.ymax - d_col)/1000;
phi_d         = phi_surf.far*exp(-phi_dcof.far*d_col);
%replicate it for entire region
phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
phi(dof_temp) = phi_total;
%Update unit classification
Mineral.unit(dof_temp) = {'Farewell'};



%Build Kma
%Locate it on the prescribed grid (get its dofs)
dof_temp      = Grid.dof(...
                (Xc > Emp_coords.ii.Kma(1) & Xc < Emp_coords.ii.Kma(2)) & ...
                (Yc > Emp_coords.ii.Kma(3) & Yc < Emp_coords.ii.Kma(4)));
dof_tot       = length(dof_temp);
d_col         = Grid.yc(Grid.yc > Emp_coords.ii.Kma(3) & ...
                        Grid.yc < Emp_coords.ii.Kma(4));
%calculate depth-dependent poro for a column of Kma
d_col         = (Grid.ymax - d_col)/1000;
phi_d         = phi_surf.ig*exp(-phi_dcof.ig*d_col);
%replicate it for entire region
phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
phi(dof_temp) = phi_total;
%Update unit classification
Mineral.unit(dof_temp) = {'Kma'};



%Build Kmrt
dof_temp      = Grid.dof(...
                (Xc > Emp_coords.ii.Kmrt(1) & Xc < Emp_coords.ii.Kmrt(2)) & ...
                (Yc > Emp_coords.ii.Kmrt(3) & Yc < Emp_coords.ii.Kmrt(4)));
dof_tot       = length(dof_temp);
d_col         = Grid.yc(Grid.yc > Emp_coords.ii.Kmrt(3) & ...
                        Grid.yc < Emp_coords.ii.Kmrt(4));
%calculate depth-dependent poro for a column of Kmrt
d_col         = (Grid.ymax - d_col)/1000;
phi_d         = phi_surf.ig*exp(-phi_dcof.ig*d_col);
%replicate it for entire region
phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
phi(dof_temp) = phi_total;
%Update unit classification
Mineral.unit(dof_temp) = {'Kmrt'};



%% 4. Edit region iii

%Build L and R
regions     = {'L','R'};
unit_names  = {'Farewell','m1','Trma','cs1','Jmdt','cs2','m2','Kmrt'};

%compaction corrected porosities for units 1-7
phi_units.L = [0.0324, 0.0406, 0.0627, 0.0467,...
               0.0732, 0.0521, 0.0574];
phi_units.R = [0.0370, 0.0451, 0.0689, 0.0518,...
               0.0804, 0.0599, 0.0638];

for i = 1:length(regions)           
    for j = 1:length(unit_names)

        %Locate it on the prescribed grid (get its dofs)
        dof_temp = Grid.dof(...
                   (Xc > Emp_coords.iii.(regions{i}).(unit_names{j})(1) & ... 
                    Xc < Emp_coords.iii.(regions{i}).(unit_names{j})(2)) & ...
                   (Yc > Emp_coords.iii.(regions{i}).(unit_names{j})(3) & ...
                    Yc < Emp_coords.iii.(regions{i}).(unit_names{j})(4)));
                
        if j == 1
            dof_tot       = length(dof_temp);
            d_col         = Grid.yc(Grid.yc > Emp_coords.iii.(regions{i}).(unit_names{j})(3) & ...
                                    Grid.yc < Emp_coords.iii.(regions{i}).(unit_names{j})(4));
            %calculate depth-dependent poro for a column of Farewell
            d_col         = (Grid.ymax - d_col)/1000;
            phi_d         = phi_surf.far*exp(-phi_dcof.far*d_col);
            %replicate it for entire region
            phi_total     = reshape(repmat(phi_d,1,dof_tot/length(d_col)),dof_tot,1);
            phi(dof_temp) = phi_total;
            %Update unit classification
            Mineral.unit(dof_temp) = unit_names(j);
        else
            %porosity is a constant value for entire region
            phi(dof_temp) = phi_units.(regions{i})(j-1);
            %Update unit classification
            Mineral.unit(dof_temp) = unit_names(j);
        end

    end
end



%Build pluton
%Locate it on the prescribed grid (get its dofs)
dof_temp = Grid.dof(...
           (Xc > Emp_coords.iii.P(1) & ... 
            Xc < Emp_coords.iii.P(2)) & ...
           (Yc > Emp_coords.iii.P(3) & ...
            Yc < Emp_coords.iii.P(4)));
%porosity is a constant value for entire pluton
phi(dof_temp) = phi_p;
%Update unit classification
Mineral.unit(dof_temp) = {'pluton'};

%Make Pluton variable (STORE DOFs)
Pluton.dof    = dof_temp;
Pluton.dof_f  = intersect(Grid.dof_f_ymin,get_dof_faces2D(Pluton.dof,Grid));

%% 5. Input vertical faults between regions

fault_thicc = 20; %[m] half thickness of fault interface

fault_bound = [Emp_coords.I(2) - fault_thicc, ...
               Emp_coords.I(2) + fault_thicc, ...
               Emp_coords.II(2) - fault_thicc, ...
               Emp_coords.II(2) + fault_thicc];
               %[m] x-coordinate limits for the fault

dof_fault   = Grid.dof((Xc > fault_bound(1) & Xc < fault_bound(2)) | ...
                       (Xc > fault_bound(3) & Xc < fault_bound(4)));
%porosity is a constant value for fault
phi(dof_fault) = phi_f;

%% Generate permeability matrix

k = perm_from_poro(k_0,phi,phi_min,phi_0,n);
k(dof_fault) = perm_from_poro(k_0,phi(dof_fault),phi_min,phi_0,n_fault);
k = reshape(k,Grid.Ny,Grid.Nx);
Kd = comp_mean(k,-1,1,Grid);

%% Create Phi field for Mineral structure based on permeability matrix

%new field to include a number identifier alongside the string
Mineral.number = zeros(Grid.N,1);

for j = 1:Grid.N
    if strcmp(Mineral.unit(j),'pluton')
        Mineral.cal.phi(j) = 0;
        Mineral.dol.phi(j) = 0;
        Mineral.qtz.phi(j) = .20;
        Mineral.An40.phi(j) = .58;
        Mineral.mus.phi(j) = 0;
        Mineral.bte.phi(j) = .15;
        Mineral.hbd.phi(j) = .07;
        Mineral.number(j)  = 1;
    elseif strcmp(Mineral.unit(j),'Farewell')
        %Arithmetic mean of phyllite and siltstone comps
        %1st term: phyllite %2nd term: siltstone
        Mineral.cal.phi(j) = .5*(0 + 0);
        Mineral.dol.phi(j) = .5*(0 + 0);
        Mineral.qtz.phi(j) = .5*(.375 + .45);
        Mineral.An40.phi(j) = .5*(.375 + .45);
        Mineral.mus.phi(j) = .5*(.25 + .10);
        Mineral.bte.phi(j) = .5*(0 + 0);
        Mineral.hbd.phi(j) = .5*(0 + 0);
        Mineral.number(j)  = 2;
    elseif any(strcmp(Mineral.unit(j),{'cs1','cs2'}))
        Mineral.cal.phi(j) = .375;
        Mineral.dol.phi(j) = .125;
        Mineral.qtz.phi(j) = .40;
        Mineral.An40.phi(j) = 0;
        Mineral.mus.phi(j) = .10;
        Mineral.bte.phi(j) = 0;
        Mineral.hbd.phi(j) = 0;
        Mineral.number(j)  = 3;
    elseif any(strcmp(Mineral.unit(j),{'m','m1','m2'}))
        Mineral.cal.phi(j) = .75;
        Mineral.dol.phi(j) = .25;
        Mineral.qtz.phi(j) = 0;
        Mineral.An40.phi(j) = 0;
        Mineral.mus.phi(j) = 0;
        Mineral.bte.phi(j) = 0;
        Mineral.hbd.phi(j) = 0;
        Mineral.number(j)  = 4;
    elseif strcmp(Mineral.unit(j),'Kmrt')
        Mineral.cal.phi(j) = 0;
        Mineral.dol.phi(j) = 0;
        Mineral.qtz.phi(j) = .55;
        Mineral.An40.phi(j) = .10;
        Mineral.mus.phi(j) = .10;
        Mineral.bte.phi(j) = .25;
        Mineral.hbd.phi(j) = 0;
        Mineral.number(j)  = 5;
    elseif any(strcmp(Mineral.unit(j),{'Jmdt','Kma'}))
        Mineral.cal.phi(j) = 0;
        Mineral.dol.phi(j) = 0;
        Mineral.qtz.phi(j) = .55;
        Mineral.An40.phi(j) = .20;
        Mineral.mus.phi(j) = .05;
        Mineral.bte.phi(j) = .15;
        Mineral.hbd.phi(j) = .05;
        Mineral.number(j)  = 6;
    end

end
Mineral.number(dof_fault) = 7;

%% Get specific points of study around Pluton

%Left side of pluton
dist_away = 50;  %[m] distance away from pluton interface where skarn forms
left  = Grid.dof(...
                 (Xc < Emp_coords.iii.P(1) & Xc > (Emp_coords.iii.P(1) - dist_away)) & ...
                  reshape((strcmp(Mineral.unit,'m2') | strcmp(Mineral.unit,'cs2')),Grid.Ny,Grid.Nx));
left  = median(left) + [1; -1];

right = Grid.dof(...
                (Xc > Emp_coords.iii.P(2) & Xc < (Emp_coords.iii.P(2) + dist_away)) & ...
                 reshape((strcmp(Mineral.unit,'m2') | strcmp(Mineral.unit,'cs2')),Grid.Ny,Grid.Nx));
right = max(right) - [1; 3];       
top   = Grid.dof(...
                (Yc > Emp_coords.iii.P(4) & ...
                 Yc < (Emp_coords.iii.P(4) + dist_away)) & ...
                 (Xc > Emp_coords.iii.P(1) & Xc < Emp_coords.iii.P(2)) & ...
                 reshape(strcmp(Mineral.unit,'m'),Grid.Ny,Grid.Nx));
top   = top(round(length(top)*[.333 .667]))';

Pluton.dof_rock = sort([left; right; top]);

end