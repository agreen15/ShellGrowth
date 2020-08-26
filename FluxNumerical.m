function [Q,z2]=FluxNumerical(T,z)   
Q=diff(T)./diff(z);
z2=(z(1:end-1)+z(2:end))/2;