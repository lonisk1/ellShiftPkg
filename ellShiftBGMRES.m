function [X,iter,rrnorm] = ellShiftBGMRES(A,B,ell,maxIter,noiseLevel,eta)

% ellShiftBGMRES - algorithm for solving block linear discrete
% problems with a square nonsymmetric matrix and initial iterate
% X_0 = [0 ... 0]. This version allows the user to select the level
% of range restriction {1,2,3,...,}. The algorithm terminates
% according to the discrepancy principle.
%
% [X,iter,rrnorm] = ellShiftBGMRES(A,B,ell,maxIter,noiseLevel,eta)
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
% GMRES, ellShiftGMRES, BGMRES, glGMRES, ellShiftGlGMRES

% Alessandro Buccini, University of Cagliari
% Lucas Onisk, Kent State University
% Lothar Reichel, Kent State University
% Code Version 1.0 - November, 2022.

breakout = eta*noiseLevel;
normB = norm(B,'fro');
[n,k] = size(B);

rrnorm = zeros(maxIter,1); %preallocate relative res. vec. 
[Q_hat,R_dag] = qr(B,0); %gives reduced QR factorization of B
H = zeros(k*(maxIter+1),k*maxIter); % preallocate block Hessenberg matrix
V = zeros(n,k*(maxIter+1)); % preallocate V matrix
V(:,1:k) = Q_hat; % assign first unitary part of V

for i = 1:ell
    W = A*V(:,k*(i-1)+1:i*k); % select proper V_i block each iteration
    for j = 1:i
        H(k*(j-1)+1:j*k,k*(i-1)+1:i*k) = V(:,k*(j-1)+1:j*k)'*W;
        W = W - V(:,k*(j-1)+1:j*k)*H(k*(j-1)+1:j*k,k*(i-1)+1:i*k);
    end
    [Q_hat,R_hat] = qr(W,0);
    V(:,i*k+1:(i+1)*k) = Q_hat; % insert new sub V_{i+1} block to V
    H(j*k+1:k*(j+1),k*(i-1)+1:i*k) = R_hat; % insert new H_{i+1,i} block to H
end

for p = 1:maxIter
    % Do next Block Arnoldi Iteration
    W = A*V(:,k*(p+ell-1)+1:(p+ell)*k);
    for j = 1:ell+p
        H(k*(j-1)+1:j*k,k*(p+ell-1)+1:(p+ell)*k) = V(:,k*(j-1)+1:j*k)'*W;
        W = W - V(:,k*(j-1)+1:j*k)*H(k*(j-1)+1:j*k,k*(p+ell-1)+1:(p+ell)*k);
    end
    [Q_hat,R_hat] = qr(W,0);
    V(:,(p+ell)*k+1:(p+ell+1)*k) = Q_hat;
    H(j*k+1:k*(j+1),k*(p+ell-1)+1:(p+ell)*k) = R_hat;
    
    % Compute necessary QR factorizations
    [Q,~] = qr(H(1:(p+1)*k,1:p*k));
    for q = 1:ell
        Q_last = Q;
        [Q,R] = qr(H(1:(q+p+1)*k,1:(q+p)*k)*Q_last(:,1:p*k));
    end
    
    % BGMRES Step
    E1 = eye((ell+p+1)*k,k);
    U = Q'*E1*R_dag;
    Y = R \ U;
    
    % Tracking
    rrnorm(p) = norm(U(p*k+1:(ell+p+1)*k,:),'fro')/normB;
    
    if rrnorm(p) <= breakout
        fprintf('**ellShiftBGMRES terminated at iteration %d** \n',p);
        X = V(:,1:(ell+p)*k) * Q_last(:,1:p*k) * Y;
        rrnorm = rrnorm(1:p);
        iter = p;
        break
    end
    if p == maxIter
        fprintf('**ellShiftedBGMRES failed to terminate after maxIters**')
        X = V(:,1:(ell+p)*k) * Q_last(:,1:p*k) * Y;
        iter = p;
        break
    end
end