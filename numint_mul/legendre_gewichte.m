function [alpha] = legendre_gewichte(n)
%LEGENDRE_GEWICHTE evaluates the legendre weigths for gauss integration of degree n

alpha=zeros(n,1);
c = legendre_roots(n);% coefficent vector
ci = zeros(n-1,1);
faktor = 0;
integral = zeros(n+1,0);
%numeric integration for polynomials
for i=1:n
    ci = [c(1:i-1)',c(i+1:n)'];
    faktor = prod(c(i)-ci);
    p = poly(ci)/faktor;%grad n-1; länge n
    integral = [p,0]./[(n:-1:1),1];
    alpha(i) = sum(integral-integral.*((-ones(1,n+1)).^(n:-1:0)));
end

end

