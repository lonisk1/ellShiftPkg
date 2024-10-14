%   blockLinearSys_demo.m (script)

%   Description: A demo to showcase the use of BGMRES, ellShiftBGMRES,
%                glGMRES, and ellShiftGlGMRES algorithms on
%                a block linear discrete ill-posed problem.

%   Instructions: Confirm that the following functions are in the working
%                 directory before running the script:
%                       - BGMRES.m
%                       - ellShiftBGMRES.m
%                       - glGMRES.m
%                       - ellShiftGlGMRES.m
%                       - phillips_alt.m

%   Functions utilized in this script:
%       [X,Iter,rrnorm] = BGMRES(A,B,maxIter,noiseLevel,eta);
%       [X,Iter,rrnorm] = ellShiftBGMRES(A,B,ell,maxIter,noiseLevel,eta);
%       [X,Iter,rrnorm] = glGMRES(A,B,maxIter,noiseLevel,eta);
%       [X,Iter,rrnorm] = ellShiftGlGMRES(A,B,ell,maxIter,noiseLevel,eta);
%       [A,b,x_true] = phillips_alt(1000);

%   Expected Results of Successful run:
%           - 4x3 table printed to command window that communicates RRE and
%             number of matrix-vector products for each method considered
%           - Figure 1 is a 2x2 plot displaying the approximate and true
%             solutions for the block methods considered for the
%             alt_phillips problem 
%%
% Clear command and workspace
clear
clc

%  Alternate_Phillips Problem by Neuman et. al
%  pkg "na33" available at http://www.netlbib.org/numeralgo/
[A,b,x_true] = phillips_alt(1000);
 
% Adjustable inputs for demo script
noiseLevel = 0.01; % 0.01 corresponds to 1% std. normal noise addition
eta = 1.01;
maxIter = 30;
ell_1 = 1; % ell value for ellShiftBGMRES function below
ell_2 = 1; % ell value for ellShiftGlobalGMRES function below

% Building block problem
r = randn(size(b));
r1 = randn(size(b));
noiseVector = ((noiseLevel*norm(b))/norm(r))*r;
noiseVector1 = ((noiseLevel*norm(b))/norm(r1))*r1;
B = [b+noiseVector b+noiseVector1];
X_true = [x_true x_true];
Xnorm = norm(X_true,'fro');

% Solution Methods
[X,iter,rrnorm] = BGMRES(A,B,maxIter,noiseLevel,1.01);
[X1,iter1,rrnorm1] = ellShiftBGMRES(A,B,ell_1,maxIter,noiseLevel,eta);
[X2,iter2,rrnorm2] = glGMRES(A,B,maxIter,noiseLevel,1.01);
[X3,iter3,rrnorm3] = ellShiftGlGMRES(A,B,ell_2,maxIter,noiseLevel,eta);

%% Table for Final RRE Values
Method = {'0-shifted BGMRES';[num2str(ell_1,'%d') '-shifted BGMRES'];'0-shifted GlGMRES';[num2str(ell_2,'%d') '-shifted GlGMRES']};
Final_RRE = [norm(X-X_true,'fro')/Xnorm;
    norm(X1-X_true,'fro')/Xnorm;
    norm(X2-X_true,'fro')/Xnorm;
    norm(X3-X_true,'fro')/Xnorm];
Iterations = [2*iter; 2*(iter1+ell_1); 2*(iter2); 2*(iter3+ell_2)]; 
        %Corresponds to number of matrix-vector products used by method

%Print table to Command Window
T = table(Method,Iterations,Final_RRE)


%% Plots

% Set-up x-axis
axis = [ones(size(A,1),1) 2.*ones(size(A,1),1)];

% Plot BGMRES Variants
figure;
subplot(2,2,1)
plot3(axis,1:size(A,1),X,'-m',axis,1:size(A,1),X_true,'-k','linewidth',1.25)
grid on
title('BGMRES soln')

subplot(2,2,2)
plot3(axis,1:size(A,1),X1,'-m',axis,1:size(A,1),X_true,'-k','linewidth',1.25)
grid on
title([num2str(ell_1,'%d') '-shifted BGMRES'])

subplot(2,2,3)
plot3(axis,1:size(A,1),X2,'-m',axis,1:size(A,1),X_true,'-k','linewidth',1.25)
grid on
title('glGMRES soln')

subplot(2,2,4)
plot3(axis,1:size(A,1),X3,'-m',axis,1:size(A,1),X_true,'-k','linewidth',1.25)
grid on
title([num2str(ell_2,'%d') '-shifted GlGMRES'])

% Figure properties
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
set(fig, 'Position',  [200, 150, 800, 600])
% Add Legend
Lgnd = legend('Block method soln','','True soln','');
Lgnd.Position(1) = 0.05;
Lgnd.Position(2) = 0.45;

