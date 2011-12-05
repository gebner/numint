function [integral, failure ] = gauss(f,a,b, epsilon )
%INTEGRIEREN_GAUSS Summary of this function goes here
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

iMin = 8;
iMax = 10;

finished = false;
failure = 1;
integral= 0;
n = iMin;

while (~finished)
    
    wert_alt = integral;
    integral = gauss_quadratur(f,n,a,b);
    failure = abs((integral-wert_alt)/integral);
    
   if (n < iMin)
       finished = 0;
   elseif(failure < epsilon)
       finished = 1;
   elseif (n > iMax)
       finished = 1;
   end
   n = n+1;
end

end

function [ wert ] = gauss_quadratur( f, n , a, b)
%GAUSS_QUADRATUR evaluates the integral using gauss integration of degree n
% @param[in] f          function pointer
% @param[in] n          degree
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[out] wert      integral value

x = legendre_roots(n);
alpha = legendre_gewichte(n);

% evaluate integrals
wert = sum((b-a)/2 *alpha.*f((b-a)/2*x+(a+b)/2));

end
