function [qx,qy,qmags] = center_flux(q,Grid)

% author: Evan J. Ramos
% date:   09 February 2016

%% Description

% This function takes the fluxes on the faces of 2D grid and applies them
% to the center of cells. The goal is to eventually make a vector field of
% new fluxes in order to see how a fluids flow will evolve through time.

%%%% NOTE: Arithmetic average of the two fluxes is utilized %%%%

qx = zeros(Grid.N,1); qy = zeros(Grid.Ny,Grid.Nx); %preallocation
oldqy = reshape(q((Grid.Nfx+1):Grid.Nf),Grid.Ny+1,Grid.Nx);

for i = 1:Grid.N
    qx(i) = (q(i+Grid.Ny) + q(i))/2;
end
qx = reshape(qx,Grid.Ny,Grid.Nx);

for j = 1:Grid.Ny
    qy(j,:) = (oldqy(j,:) + q(j+1,:))/2;
end

% Magnitude
qmags = sqrt(qx.^2 + qy.^2);

end