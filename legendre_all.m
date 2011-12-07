function [ x, alpha ] = legendre_all( iMin, iMax )
%LEGENDRE_ALL derives the roots of the legendre polynomial and the weights
%   Detailed explanation goes here

x = zeros(iMax);
alpha = zeros(iMax);

for i = iMin:iMax
    x(i,(iMax-i+1):iMax) = legendre_roots(i);
    alpha(i,(iMax-i+1):iMax) = legendre_gewichte(i);
end

end

