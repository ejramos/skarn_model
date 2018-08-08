function u = solve_lbvp(L,fs,B,g,N) 


%% Q2

% author: Evan Ramos 
% date:   20 September 2015 

% Description: 
% Computes the solution $u$ to the linear differential problem given by 
% % $$\mathcal{L}(u)=f \quad x\in \Omega $$ 
% % with boundary conditions 

% % $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$. 

% % Input: 
% L = matrix representing the discretized linear operator of size N by N, 
% where N is the number of degrees of fredom 
% f = column vector representing the discretized r.h.s. and contributions 
% due non-homogeneous Neumann BC’s of size N by 1 
% B = matrix representing the constraints arising from Dirichlet BC’s of 
% size Nc by N 
% g = column vector representing the non-homogeneous Dirichlet BC’s of size 
% Nc by 1. 
% N = matrix representing a orthonormal basis for the null-space of B and 
% of size N by (N-Nc). 

% Output: 
% u = column vector of the solution of size N by 1
if isempty(B)
    u = L\fs;
else
    hp = B'*((B*B')^-1)*g; %particular solution
    u = N*((N'*L*N)\(N'*(fs-L*hp)))+hp; %homogeneous solution
end


end