function [ x, alpha ] = gauss_arrays( iMin, iMax )
%LEGENDRE_ALL derives the roots of the legendre polynomial and the weights
% @param[in]    iMin        minimal degree of the legendre polynomial
% @param[in]    iMax        maximal degree of the legendre polynomial
% @param[out]   x           matrix of the legendre roots (for degree i last
% i elements of the i-th row)
% @param[out]   alpha       matrix of the weights

x = zeros(iMax);
alpha = zeros(iMax);

for i = iMin:iMax
    x(i,(iMax-i+1):iMax) = legendre_roots(i);
    alpha(i,(iMax-i+1):iMax) = legendre_gewichte(i);  
end

end