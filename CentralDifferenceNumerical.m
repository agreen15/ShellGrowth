function [Q,z2]=CentralDifference(T,z)   
Q=(T(3:end)-T(1:end-2))./(z(3:end)-z(1:end-2));
z2=(z(3:end)+z(1:end-2))/2;