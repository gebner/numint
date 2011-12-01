function [ wert ] = gauss_quadratur( f, n , a, b)
%GAUSS_QUADRATUR Summary of this function goes here
%   Detailed explanation goes here
x = legendre_roots(n);
alpha = legendre_gewichte(n);

% wert = 0;
% for i = 1:n
%     wert = wert + alpha(i)*f((b-a)/2*x(i)+(a+b)/2);
% end
% wert = wert*(b-a)/2;
wert = sum((b-a)/2 *alpha.*f((b-a)/2*x+(a+b)/2));
end

