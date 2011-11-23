function[I] = intpl_int(f, a, b, n, intpl_worker)
% intpl_int Evaluate the integral of the intepolation-polynomial of f
% returned by intpl_worker on the interval [a, b].
%
% Example call:
% intpl_int(@(x)exp(x), 0, 1, 3, @(f, a, b, n)aqvdst_pf(f, a, b, n))
%
% f ... function
% a ... lower interval bound
% b ... upper interval bound
% n ... degree of the interpolation
% intpl_worker ... an interpolation function that takes f, a, b and n
% as arguments.

    poly = intpl_worker(f, a, b, n);
    poly = polyint(poly);
    I = polyval(poly, b) - polyval(poly, a);
end
