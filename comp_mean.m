function [Kd] = comp_mean(K,p,kvkh,Grid) 

%% Q2: Part 1

% author: Evan Ramos
% date:   14 October 2015

% Description: 
% Takes coefficient field, K, defined at the cell centers and computes the 
% mean specified by the power, p and returns it in a sparse diagonal 
% matrix, Kd. 


% % Input: 
% K = Ny by Nx matrix of cell centered values 
% p = power of the generalized mean 
%          1 (arithmetic mean) 
%         -1 (harmonic mean)
% kvkh = ratio of vertical to horizontal conductivity/permeability
%        (anisotropy)
% Grid = structure containing information about the grid. 


% % Output: 
% Kd = Nf by Nf diagonal matrix of power means at the cell faces. 


% % Example call: 
% K = @(x) 1+x.^3; 
% Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% Grid = build_grid(Grid); 
% Kd = comp_mean(K(Grid.xc),1,Grid);


if (p == -1) || (p == 1)
    
    if (Grid.Nx == Grid.N) || (Grid.Ny == Grid.N) % 1D
        mean = zeros(Grid.Nx+1,1);
        mean(2:Grid.Nx) = sum(.5*[K(1:Grid.Nx-1),K(2:Grid.Nx)].^p,2).^(1/p);
        Kd = spdiags(mean,0,Grid.Nx+1,Grid.Nx+1);
        
    elseif (Grid.N > Grid.Nx) || (Grid.N > Grid.Ny) % 2D
        mean_x = zeros(Grid.Ny,Grid.Nx+1);
        mean_x(:,2:Grid.Nx) = ...
            (.5*(K(:,1:Grid.Nx-1)+K(:,2:Grid.Nx)).^p).^(1/p);
        mean_y = zeros(Grid.Ny+1,Grid.Nx);
        mean_y(2:Grid.Ny,:) = ...
            (.5*(K(1:Grid.Ny-1,:)+K(2:Grid.Ny,:)).^p).^(1/p);
        onevec = [reshape(kvkh*mean_x,1,numel(mean_x)),...
            reshape((1/kvkh)*mean_y,1,numel(mean_y))];
        Kd = spdiags(onevec',0,Grid.Nf,Grid.Nf);
    else
        error('3D permeability is not implemented')
    end
    
else
    error('This power does not have significance.')
end

end