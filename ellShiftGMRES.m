function [x,iter,rrnorm] = ellShiftGMRES(A,b,ell,maxIter,noiseLevel,eta)

% ellShiftGMRES - algorithm for solving linear discrete problems with a
% square nonsymmetric matrix and initial iterate x_0 = 0. This version 
% allows the user to select the level of range restriction {1,2,3,...}. 
% The algorithm terminates according to the discrepancy principle.
%
% [x,iter,rrnorm] = ellShiftGMRES(A,b,ell,maxIter,noiseLevel,eta)
%
% Inputs: A - The (n x n) matrix of the linear discrete ill-posed problem Ax = b.
%         b - The right-hand side of the problem Ax = b.
%         ell - A positive integer that defines the shift of the Krylov
%         subspace where a solution is sought.
%         maxIter - The maximum number of iterations the algorithm
%                   will attempt to carryout if the discrepancy
%                   principle is not satisifed.
%         noiseLevel - The scaled noise level contaminating the b
%                       vector. The scaled noise level is given by
%                       100(||e||/||b||) where ||e|| is the norm of the
%                       error contaminating the right-hand side
%                       vector b. Typical values are 0.01 (i.e. 1%
%                       noise), 0.005 (i.e. 0.5% noise), etc.
%         eta - A constant > 1 used in the discrepancy principle (DP).
%               Standard usage is eta = 1.01.
%
% Outputs: x - Approximate solution to Ax = b.
%          iter - The number of computed iterations for the algorithm.
%          rrnorm - The relative residual norm vector for each
%                   iteration given by ||b - Ax_i|| / ||b||
%                   for i = 1,2,....
%
% See also:
% GMRES, BGMRES, ellShiftBGMRES, glGMRES, ellShiftGlGMRES

% Alessandro Buccini, University of Cagliari
% Lucas Onisk, Kent State University
% Lothar Reichel, Kent State University
% Code Version 1.0 - November, 2022.

breakout = eta*noiseLevel;
normB = norm(b,2);

rrnorm = zeros(maxIter,1); %preallocate relative res. vec.
V = zeros(size(b,1),ell+maxIter+1); 
V(:,1) = b/norm(b,2); %start V matrix
H = zeros(ell+maxIter,ell+maxIter+1); %preallocate Hessenberg matrix

for i = 1:ell
    v = A*V(:,i);
    for j = 1:i
        H(j,i) = V(:,j)'*v;
        v = v - H(j,i)*V(:,j);
    end
    H(i+1,i) = norm(v,2);
    V(:,i+1) = v/H(i+1,i); 
end

for p = 1:maxIter
    % Do next Arnoldi Iteration
    v = A*V(:,p+ell);
    for j = 1:p+ell
        H(j,p+ell) = V(:,j)'*v;
        v = v - H(j,p+ell)*V(:,j);
    end
    H(p+ell+1,p+ell) = norm(v,2);
    V(:,ell+p+1) = v/H(p+ell+1,p+ell); 
    
    % Compute necessary QR factorizations
    [Q,~] = qr(H(1:p+1,1:p));
    for j = 1:ell
        Q_last = Q;
        [Q,R] = qr(H(1:j+p+1,1:j+p)*Q_last(:,1:p));
    end
    
    e1 = zeros(size(Q,1),1);
    e1(1) = normB;
    g = Q'*e1;
    y = R \ g;
    
    % Tracking residual
    rrnorm(p) = norm(g(p+1:length(g)),2)/normB;
    
    if rrnorm(p) < breakout
        fprintf('**ellShiftGMRES terminated at iteration %d** \n',p);
        rrnorm = rrnorm(1:p);
        x = V(:,1:ell+p) * Q_last(:,1:p) * y; %approximate soln
        iter = p;
        break
    end
    if p == maxIter
        fprintf('**ellShiftedGMRES failed to terminate after maxIters**')
        x = V(:,1:ell+p) * Q_last(:,1:p) * y; %approximate soln
        iter = p;
        break
    end
end