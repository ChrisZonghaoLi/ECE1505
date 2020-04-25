close all
clear all
syms b
p = ((2*b-2)/(4-b^2))^2 + ((b-4)/(4-b^2))^2 + b*((2*b-2)/(4-b^2))*...
    ((b-4)/(4-b^2)) + (2*b-2)/(4-b^2) + (b-4)/(4-b^2);
expand(p);
p_prime = diff(p)

eqn = p_prime == 0;
beta = solve(eqn)