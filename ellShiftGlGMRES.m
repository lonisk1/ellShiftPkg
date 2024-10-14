function [X,iter,rrnorm] = ellShiftGlGMRES(A,B,ell,maxIter,noiseLevel,eta)

% ellShiftGlGMRES - algorithm for solving block linear discrete
% problems with a square nonsymmetric matrix and initial iterate
% X_0 = [0 ... 0]. This version allows the user to select the level
% of range restriction {1,2,3,...,}. The algorithm terminates
% according to the discrepancy principle.
%
% [X,iter,rrnorm] = ellShiftGlGMRES(A,B,ell,maxIter,noiseLevel,eta)
%
% Inputs: A - The (n x n) matrix of the linear discrete ill-posed problem AX = B.
%         B - The right-hand side (n x k) matrix of the problem AX = B.
%         ell - A positive integer that defines the shift of the Krylov
%         subspace where a solution is sought.
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
% GMRES, ellShiftGMRES, BGMRES, ellShiftBGMRES, glGMRES

% Alessandro Buccini, University of Cagliari
% Lucas Onisk, Kent State University
% Lothar Reichel, Kent State University
% Code Version 1.0 - November, 2022.

breakout = eta*noiseLevel;
normB = norm(B,'fro');
[n,k] = size(B);

rrnorm = zeros(maxIter,1); %preallocate relative res. vec. 
V = zeros(n,k*(maxIter+1)); % preallocate V matrix
V(:,1:k) = B/normB;
H = zeros((maxIter+1),maxIter); % preallocate Hessenberg matrix

for i = 1:ell
    W = A*V(:,k*(i-1)+1:i*k); % select proper V_i block each iteration
    for j = 1:i
        standIn = V(:,k*(j-1)+1:j*k);
        H(j,i) = W(:)'*standIn(:);
        W = W - H(j,i).*V(:,k*(j-1)+1:j*k);
    end
    H(j+1,i) = norm(W,'fro');
    V(:,i*k+1:(i+1)*k) = W./H(j+1,i);
end

for p = 1:maxIter
    % Do next Globl Arnoldi Iteration
    W = A*V(:,k*(p+ell-1)+1:(p+ell)*k);
    for j = 1:ell+p
        standIn = V(:,k*(j-1)+1:j*k);
        H(j,p+ell) = W(:)'*standIn(:);
        W = W - H(j,p+ell).*standIn;
    end
    
    H(j+1,p+ell) = norm(W,'fro');
    V(:,(p+ell)*k+1:(p+ell+1)*k) = W./H(j+1,p+ell);
    
    % Compute necessary QR factorizations
    [Q,~] = qr(H(1:p+1,1:p));
    for q = 1:ell
        Q_last = Q;
        [Q,R] = qr(H(1:q+p+1,1:q+p)*Q_last(:,1:p));
    end
    
    e1 = eye(ell+p+1,1);
    U = Q'*e1.*normB;
    Y = R \ U;
    
    % Tracking
    rrnorm(p) = norm(U(p+1:length(U)),2)/normB;
    
    if rrnorm(p) <= breakout
        fprintf('**ellShiftGlGMRES terminated at iteration %d** \n',p);
        X = V(:,1:(ell+p)*k) * kron((Q_last(:,1:p)*Y),eye(k)); %approx. soln
        rrnorm = rrnorm(1:p);
        iter = p;
        break
    end
    if p == maxIter
        fprintf('ellShiftedGlGMRES failed to terminate after maxIters')
        X = V(:,1:(ell+p)*k) * kron((Q_last(:,1:p)*Y),eye(k)); %approx. soln
        iter = p;
        break
    end
end