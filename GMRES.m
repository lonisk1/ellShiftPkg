function [x,iter,rrnorm] = GMRES(A,b,maxIter,noiseLevel,eta)

% GMRES - algorithm for solving linear discrete problems with a
% square nonsymmetric matrix and initial iterate x_0 = 0.
% The algorithm terminates according to the discrepancy principle.
%
% [x,iter,rrnorm] = GMRES(A,b,maxIter,noiseLevel,eta)
%
% Inputs: A - The (n x n) matrix of the linear discrete ill-posed problem Ax = b.
%         b - The right-hand side of the problem Ax = b.
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
% ellShiftGMRES, BGMRES, ellShiftBGMRES, glGMRES, ellShiftGlGMRES

% Alessandro Buccini, University of Cagliari
% Lucas Onisk, Kent State University
% Lothar Reichel, Kent State University
% Code Version 1.0 - November, 2022. 

breakout = eta*noiseLevel; % "breakout" according to DP
normB = norm(b,2);

rrnorm = zeros(maxIter,1); %preallocate relative res. vec.
Q = zeros(size(b,1),maxIter+1);
Q(:,1) = b/norm(b,2); %assign first column of Q matrix
H = zeros(maxIter,maxIter+1); %preallocate Hessenberg matrix

for i = 1:maxIter
    v = A*Q(:,i);
    for j = 1:i
        H(j,i) = Q(:,j)'*v;
        v = v - H(j,i)*Q(:,j);
    end
    
    H(i+1,i) = norm(v,2);
    Q(:,i+1) = v/H(i+1,i);
    
    % GMRES Step
    [q,r] = qr(H(1:i+1,1:i));
    g = normB.*(q'*eye(i+1,1));
    y = r\g;
    
    % Tracking residual
    rrnorm(i) = abs(g(i+1))/normB;
    
    if rrnorm(i) < breakout
        fprintf('**GMRES terminated at iteration %d** \n',i);
        rrnorm = rrnorm(1:i);
        x = Q(:,1:i)*y; %approximate soln
        iter = i;
        break
    end
    if i == maxIter
        fprintf('**GMRES failed to terminate after maxIters')
        x = Q(:,1:i)*y; %approximate soln
        iter = i;
        break
    end
end