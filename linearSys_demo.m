%   linearSys_demo.m (script)

%   Description: A demo to showcase the use of GMRES and ellShiftGMRES
%                algorithms on a linear discrete ill-posed problem.

%   Instructions: Confirm that the following functions are in the working
%                 directory before running the script:
%                       - GMRES.m
%                       - ellShiftGMRES.m
%                       - shaw_alt.m

%   Functions utilized in this script:
%       [X,Iter,rrnorm] = GMRES(A,bn,maxIter,noiseLevel,eta);
%       [X,Iter,rrnorm] = ellShiftGMRES(A,bn,ell,maxIter,noiseLevel,eta);
%       [A,b,x_true] = shaw_alt(1000);

%   Expected Results of Successful run:
%           - 4x3 table printed to command window that replicates results
%             found in manuscript
%           - Figure 1 is a 1x4 plot displaying the true solution for the 
%             alt_shaw problem compared against the computed solution for
%             the specfied method
%           - Figure 2 displays the relative residual evolution of the
%           considered methods
%%
% Clear command and workspace
clear
clc

%  Alternate_Shaw Problem by Neuman et. al
%  pkg "na33" available at http://www.netlbib.org/numeralgo/
[A,b,x_true] = shaw_alt(1000);

% Adjustable inputs for demo script
noiseLevel = 0.01; % 0.01 corresponds to 1% std. normal noise addition
eta = 1.01;
maxIter = 30;

% building linear problem
seed = rng(6,'philox'); %setting seed number for random number generator
r = randn(size(b));
noiseVector = ((noiseLevel*norm(b))/norm(r))*r;
bn = b + noiseVector;

% Solution Methods
[X,Iter,rrnorm] = GMRES(A,bn,maxIter,noiseLevel,eta);
[X1,Iter1,rrnorm1] = ellShiftGMRES(A,bn,1,maxIter,noiseLevel,eta);
[X2,Iter2,rrnorm2] = ellShiftGMRES(A,bn,2,maxIter,noiseLevel,eta);
[X3,Iter3,rrnorm3] = ellShiftGMRES(A,bn,3,maxIter,noiseLevel,eta);

%% Table for Final RRE Values
Method = {'0-shifted GMRES';'1-shifted GMRES';'2-shifted GMRES';'3-shifted GMRES'};
Final_RRE = [norm(X-x_true,2)/norm(x_true,2);
    norm(X1-x_true,2)/norm(x_true,2);
    norm(X2-x_true,2)/norm(x_true,2);
    norm(X3-x_true,2)/norm(x_true,2)];
Iterations = [Iter; Iter1+1; Iter2+2; Iter3+3]; 
        %Corresponds to number of matrix-vector products used by method

%Print table to Command Window
T = table(Method,Iterations,Final_RRE)

%% Plots
gcf = figure(1);
set(gcf, 'Position',  [40, 200, 1200, 400])
subplot(1,4,1)
plot(1:length(x_true),x_true,'-k',1:length(X),X,'-g','linewidth',0.8)
legend('True Soln','0-shift Soln')

subplot(1,4,2)
plot(1:length(x_true),x_true,'-k',1:length(X1),X1,'-m','linewidth',0.8)
legend('True Soln','1-shift Soln')

subplot(1,4,3)
plot(1:length(x_true),x_true,'-k',1:length(X2),X2,'-r','linewidth',0.8)
legend('True Soln','2-shift Soln')

subplot(1,4,4)
plot(1:length(x_true),x_true,'-k',1:length(X3),X3,'-c','linewidth',0.8)
legend('True Soln','3-shift Soln')

figure(2);
breakoutVec = (noiseLevel*eta).*ones(max([length(rrnorm) length(rrnorm1) length(rrnorm2) length(rrnorm3)]),1);
semilogy(1:length(breakoutVec),breakoutVec,'--k',1:length(rrnorm),rrnorm,'-og',1:length(rrnorm1),rrnorm1,'-^m',1:length(rrnorm2),rrnorm2,'-rx',1:length(rrnorm3),rrnorm3,'-*c','linewidth',0.8)
legend('breakout','0-shift','1-shift','2-shift','3-shift')
title('Relative Residual Plot')