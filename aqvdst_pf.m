function[pol] = aqvdst_pf(f, a, b, n)
% aqvdst_pf Interpolates the function f through a polynomial of degree n
% returned by the matlab function polyfit given equidistant evaluation
% points on the interval [a,b].
%
% f ... function
% a ... lower interval bound
% b ... upper interval bound
% n ... degree of the interpolation

    x = a:(b-a)/n:b;
    y = arrayfun(@(x) f(x), x);
    pol = polyfit(x,y,n);
end
