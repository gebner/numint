function [ c ] = legendre_polynomial( n )
%LEGENDRE_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here
c1 = legendre_pol_ganz(n);
c = c1(n+1,1:n+1);
end

