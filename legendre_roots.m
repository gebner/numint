function [ x ] = legendre_roots(n)
%LEGENDRE_ROOTS evaluates the roots of the n-th legendre polynomial

c = legendre_polynomial(n);
x = roots(c);

end

%----------------------------------------------------------------------------------
function [ c ] = legendre_polynomial(n)
%LEGENDRE_POLYNOMIAL derives the coefficents of the n-th legendre polynomial

c1 = legendre_pol_ganz(n);
c = c1(n+1,1:n+1);

end


%-----------------------------------------------------------------------
function [c ] = legendre_pol_ganz( n)
%LEGENDRE_POL_GANZ evaluates the legendre coefficents by recursion 

if (n==0)
    c=[1];
elseif (n==1)
    c=[0,1;
       1,0];
else
    c1 = legendre_pol_ganz(n-1);
    c =[zeros(n,1),c1;(2*n-1)/n*[c1(n,1:n),0]-(n-1)/n*[0,c1(n-1,1:n)]];
end

end


