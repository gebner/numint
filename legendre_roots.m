function [ x ] = legendre_roots(n)
%LEGENDRE_ROOTS Summary of this function goes here
%   Detailed explanation goes here
c = legendre_polynomial(n);

x = roots(c);

end

