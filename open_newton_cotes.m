% I = open_newton_cotes(f, a,b, degree, num_points)
%
% Integrates the function f on the interval [a,b] using the Newton-Cotes
% formula of the specified degree and (approximately) num_points evaluation
% points.
%
function I = open_newton_cotes(f, a,b, d, n),
n = ceil((n-1)/(d+2));
w = (d+2) * (vander((1:d+1)/(d+2))' \ (1./fliplr(1:d+1))')';
w = [repmat([0 w], 1, n) 0];
x = linspace(a,b, length(w));
h = x(2) - x(1);
I = h * w * f(x)';
