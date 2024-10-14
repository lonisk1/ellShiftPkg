function [X,iter,rrnorm] = BGMRES(A,B,maxIter,noiseLevel,eta)

% BGMRES - algorithm for solving block linear discrete problems
% with a square nonsymmetric matrix and initial iterate X_0 = [0 ... 0]. 
% The algorithm terminates according to the discrepancy principle.
%
% [X,iter,rrnorm] = BGMRES(A,B,maxIter,noiseLevel,eta)
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
% Outputs: X - Approximate solution to AX = B.
%          iter - The number of computed iterations for the algorithm.
%          rrnorm - The relative residual norm vector for each
%                   iteration given by ||B - AX_i|| / ||B||
%                   for i = 1,2,....
%
% See also:
% GMRES, ellShiftGMRES, ellShiftBGMRES, glGMRES, ellShiftGlGMRES

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
Q = zeros(n,k*(maxIter+1)); % preallocate Q matrix
Q(:,1:k) = Q_hat; % assign first unitary part of Q

for j = 1:maxIter
    W = A*Q(:,k*(j-1)+1:j*k); % select proper Q_j block each iteration
    for i = 1:j
        H(k*(i-1)+1:i*k,k*(j-1)+1:j*k) = Q(:,k*(i-1)+1:i*k)'*W;
        W = W - Q(:,k*(i-1)+1:i*k)*H(k*(i-1)+1:i*k,k*(j-1)+1:j*k);
    end
    [Q_hat,R_hat] = qr(W,0);
    Q(:,j*k+1:(j+1)*k) = Q_hat; % insert new sub Q_{j+1} block to Q
    H(i*k+1:k*(i+1),k*(j-1)+1:j*k) = R_hat; % insert new H_{j+1,j} block to H
    
    % Block GMRES Step
    [q,r] = qr(H(1:(j+1)*k,1:j*k));
    U = q'*eye(k*(i+1),k)*R_dag;
    Y = r\U;
    
    % Tracking
    rrnorm(j) = norm(U(j*k+1:(j+1)*k,:),'fro')/normB;
    
    if rrnorm(j) <= breakout
        fprintf('***BGMRES Terminated at iteration %d*** \n',j);
        X = Q(:,1:j*k)*Y; %approximate soln
        rrnorm = rrnorm(1:j);
        iter = j;
        break
    end
    if j == maxIter
        fprintf('**BGMRES failed to terminate after maxIters**')
        X = Q(:,1:j*k)*Y; %approximate soln
        iter = j;
        break
    end
end