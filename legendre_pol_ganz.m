function [c ] = legendre_pol_ganz( n)
%LEGENDRE_POL_GANZ Summary of this function goes here
%   Detailed explanation goes here
if(n==0)
    c=[1];
elseif(n==1)
    c=[0,1;
       1,0];
else
    %c=(2*n-1)/(n)*[legendre_polynomial(n-1),0]-(n-1)/n*[0,0,legendre_polynomial(n-2)];
    c1 = legendre_pol_ganz(n-1);
    %c2 = [zeros(n,1),c1];
    c =[zeros(n,1),c1;(2*n-1)/n*[c1(n,1:n),0]-(n-1)/n*[0,c1(n-1,1:n)]];
end

end

