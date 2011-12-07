function [integral, failure] = gauss(f,a,b, epsilon, x, alpha, iMin, iMax)
%INTEGRIEREN_GAUSS Summary of this function goes here
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[in] x          matrix of legendre roots (optional)
% @param[in] alpha      matrix of legendre weights (optional)
% @param[in] iMin       minimum index (optional)
% @param[in] iMax       maximum index (optional)
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

if(nargin < 8)
    iMin = 8;
    iMax = 10;
end
if (nargin < 6)
    [x, alpha] = gauss_arrays(iMin, iMax);
end

failure = Inf;
integral= Inf;
i = iMin;

while (failure > epsilon && i <= iMax)
    
   wert_alt = integral;
   %evaluate integral
   integral = gauss_quadratur(f,a,b,x,alpha,i);
   % get failure estimation
   failure = abs((integral - wert_alt)/integral);
   
   i = i+1;
end

end

function [ wert ] = gauss_quadratur( f, a, b, x, alpha, n)
%GAUSS_QUADRATUR evaluates the integral using gauss integration of degree n
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] x          matrix of legendre roots
% @param[in] alpha      matrix of legendre weights
% @param[in] n          degree
% @param[out] wert      integral value

% evaluate integrals
wert = sum((b-a)/2 *[alpha(n,(end-n+1):end)].*f((b-a)/2*[x(n,(end-n+1):end)]+(a+b)/2));

end
