function [integral, failure, A] = gauss(f,a,b, epsilon, A, iMin, iMax, x, alpha)
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
% @param[out] A         statistics array
if (nargin <5)
    A=zeros(2,4);
    f = @(x) mapWithStats(f,x); % add dummy stats array
end

if(nargin < 7)
    iMin = 8;
    iMax = 10;
end
if (nargin < 9)
    [x, alpha] = gauss_arrays(iMin, iMax+1);
end
failure = Inf;
integral= Inf;
i = iMin;

while (failure > epsilon && i <= iMax)
   wert_alt = integral;
   %evaluate integral
   [integral,Err,At] = gauss_quadratur(f,a,b,x,alpha,i);
   A = A+At; %refresh statistics
   
   % get failure estimation
   if(integral~=0)
        failure = abs((abs(integral - wert_alt)+Err)/integral);
   else
       failure = abs(abs(integral-wert_alt)+Err);
   end
   
   i = i+1;
end

end

function [ I , Err,A] = gauss_quadratur( f, a, b, x, alpha, n)
%GAUSS_QUADRATUR evaluates the integral using gauss integration of degree n
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] x          matrix of legendre roots
% @param[in] alpha      matrix of legendre weights
% @param[in] n          degree
% @param[out] I         integral value
% @param[out] Err       error approximation
% @param[out] A         statistics array

[I,A]=f((b-a)/2*x(n,(end-n+1):end)+(a+b)/2);
s=size(I);
if s(1)==2
    Err=abs(I(2,1:end).*I(1,1:end));
    Err=abs((b-a)/2 *alpha(n,(end-n+1):end))*Err';
else
    Err=0;
end
I = (b-a)/2 *alpha(n,(end-n+1):end) * I(1,1:end)';
end
