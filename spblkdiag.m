function [Ablk,An] = spblkdiag(A,N)

%% Q2

% author: Evan Ramos
% date:   18 October 2015

% Description:
% This function builds a sparse block diagonal matrix, Ablk, with N
% identical blocks A.

% Input:
% A = sparse matrix.
% N = number of A-blocks in the sparse block diagonal matrix Ablk.

% Output:
% Ablk = sparse block diagonal matrix with N identical 
%        blocks corresponding to A

% Example call:
% >> A = [1 0 2; 0 0 3]; A = sparse(A);
% >> Ablk = spblkdiag(A,5)

if ~iscell(A)
    An = repmat(A,1,N);     % replicating A, N times
    [r,c] = size(A);        % size of original sparse matrix A

    An = mat2cell(An,r,c*ones(1,N));
else
    An = A;
end

Ablk = blkdiag(An{:});
end