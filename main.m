%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Bézier walk evolution (BWE) source codes version 1.0
%  
%  Developed in:	MATLAB 23.2.0.2365128 (R2023b)
%  
%  Programmer:		Jinpeng Wang
%                   e-mail: wangjinpengchunuo@163.com &
%                   2211060117@stu.lntu.edu.cn
%  
%  Original paper:	Jinpeng Wang, Xingguo Xu, Yujing Sun, Jiguang Yu,
%                   Kaichen Ouyang, and Yuansheng Gao.
%                   Random Walk on Bézier Curves for
%                   Global Optimization
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To use this code in your own project 
% remove the line for 'CEC2017' function 
% and define the following parameters: 
% fun   : function handle to the .m file containing the objective function
%		  the .m file you define should accept the whole population 'x' 
%		  as input and return a column vector containing objective function 
%		  values of all of the population members
% nvars : number of decision/design variables 
% lb    : lower bound of decision variables (must be of size 1 x nvars)
% ub    : upper bound of decision variables (must be of size 1 x nvars)
%
% BWE will return the following: 
% x     : best solution found 
% fval  : objective function value of the found solution 

clc
clear
close all

%% Inputs 
FunctionNumber = 'F1'; % F1~F30 except F2
D = 50;
N = 50;
T  = 500;
[lb,ub,nvars,fun] = GetFunctionsDetails(FunctionNumber, D);

%% Run BWE
[x,fval,ConvergenceCurve] = BWE (fun, nvars, lb, ub, N, T);
%% Plot results 
figure
plot(log(ConvergenceCurve));
xlabel('Iterations');
ylabel('Ln Objective function');
title('Convergence curve');