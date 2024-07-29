clc 
close all
clear 

%% Optimization using Genetic Algorithm (GA)

x_rec = zeros(5,5) ; 

for m = 1:5 % Perform the GA 5 times

nvars = 5; % Number of variables
A = [];
B = [];
Aeq = []; beq = []; intcon = []; nonlcon = [];
lb = [0.66, 0.92, 88 , 88 , 88 ]; % lower bounds of a, b, aplph, zeta and gamma
ub = [0.72, 0.98, 92, 92 , 92 ]; % upper bounds of a, b, aplph, zeta and gamma

options = optimoptions('ga', ...                                    
                       'PopulationSize', 10, ...                    
                       'MaxGenerations', 1000, ...                   
                       'CrossoverFraction', 0.7, ...                
                       'MaxStallGenerations', 20,...                           
                       'Display', 'iter', ...                          
                       'PlotFcn', {@gaplotbestf});    
        

[x,fval] = ga(@(x) Kossels(x), nvars, A, B, Aeq, beq, lb, ub, nonlcon, intcon, options);

x_rec(m,:) = x ;

end