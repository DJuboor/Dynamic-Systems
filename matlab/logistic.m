function [ f ] = logistic( t,t0,L,k )
% The logistic function.
% L is the limiting maximum value of the sigmoidal curve,
% t0 marks the midpoint and k is the steepness. 

f = L ./ (1 + exp(-k*(t-t0)));

end

