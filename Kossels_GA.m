clc 
close all
clear 

x_rec = zeros(5,5) ; 

for m = 1:5

%% Optimization using Genetic Algorithm (GA)
nvars = 5;
A = [];
B = [];
Aeq = []; beq = []; intcon = []; nonlcon = [];
lb = [0.66, 0.92, 88 , 88 , 88 ]; % lower bound: x = [a b strain_angle]
ub = [0.72, 0.98, 92, 92 , 92 ]; % upper bound: x = [a b strain_angle]

options = optimoptions('ga', ...                                    
                       'PopulationSize', 10, ...                    
                       'MaxGenerations', 1000, ...                   
                       'CrossoverFraction', 0.7, ...                
                       'MaxStallGenerations', 20,...                           
                       'Display', 'iter', ...                          
                       'PlotFcn', {@gaplotbestf});    
        

[x,fval] = ga(@(x) Kossels(x), nvars, A, B, Aeq, beq,  lb, ub, nonlcon, intcon, options);
% [x,fval] = ga(@(x) BTN(x,H1,Dn1,IP,lambda,LimitPrecision),nvars,[],[],[],[],lb,ub)

x_rec(m,:) = x ;

end