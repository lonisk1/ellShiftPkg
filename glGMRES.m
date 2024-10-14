function [X,iter,rrnorm] = glGMRES(A,B,maxIter,noiseLevel,eta)

% glGMRES - algorithm for solving block linear discrete
% problems with a square nonsymmetric matrix and initial 
% iterate X_0 = [0 ... 0]. The algorithm terminates
% according to the discrepancy principle.
%
% [X,iter,rrnorm] = glGMRES(A,B,maxIter,noiseLevel,eta)
%
% Inputs: A - The (n x n) matrix of the linear discrete ill-posed problem AX = B.
%         B - The right-hand side (n x k) matrix of the problem AX = B.
%         maxIter - The maximum number of iterations the algorithm
%                   will attempt to carryout if the discrepancy
%                   principle is not satisifed. 
%         noiseLevel - The scaled noise level contaminating the
%                       right-hand side B. The scaled noise level is given by
%                       100(||e||/||B||) where ||e|| is the norm of the
%                       error contaminating the right-hand side
%                       matrix B. Typical values are 0.01 (i.e. 1%
%                       noise), 0.005 (i.e. 0.5% noise), etc.
%         eta - A constant > 1 used in the discrepancy principle (DP).
%               Standard usage is eta = 1.01.
%
% Outputs: x - Approximate solution to AX = B.
%          iter - The number of computed iterations for the algorithm.
%          rrnorm - The relative residual norm vector for each
%                   iteration given by ||B - AX_i|| / ||B||
%                   for i = 1,2,....
%
% See also:
% GMRES, ellShiftGMRES, BGMRES ,ellShiftBGMRES, ellShiftGlGMRES

% Alessandro Buccini, University of Cagliari
% Lucas Onisk, Kent State University
% Lothar Reichel, Kent State University
% Code Version 1.0 - November, 2022.

normB = norm(B,'fro');
breakout = eta*noiseLevel;
[n,k] = size(B);

rrnorm = zeros(maxIter,1); %preallocate relative res. vec. 
V = zeros(n,k*(maxIter+1)); % preallocate V matrix
V(:,1:k) = B/normB; % assign first unitary part of V
H = zeros(maxIter+1,maxIter); % preallocate Hessenberg matrix

for j = 1:maxIter
    W = A*V(:,k*(j-1)+1:j*k);
    for i = 1:j
        standIn = V(:,k*(i-1)+1:i*k);
        H(i,j) = W(:)'*standIn(:);
        W = W - H(i,j).*V(:,k*(i-1)+1:i*k);
    end
    H(i+1,j) = norm(W,'fro');
    V(:,j*k+1:k*(j+1)) = W./H(i+1,j);
    
    % Global GMRES Step
    [Q,R] = qr(H(1:j+1,1:j));
    g = normB.*(Q'*eye(j+1,1));
    y = R\g;
    
    rrnorm(j) = abs(g(j+1))/normB;
    if rrnorm(j) <= breakout
        fprintf('**glGMRES terminated at iteration %d** \n',j);
        X = V(:,1:j*k)*kron(y,eye(k)); %approx. soln.
        rrnorm = rrnorm(1:j);
        iter = j;
        break;
    end
    if j == maxIter
        fprintf('**glGMRES failed to terminate after maxIters**')
        X = V(:,1:j*k)*kron(y,eye(k)); %approx. soln.
        iter = j;
        break
    end
end