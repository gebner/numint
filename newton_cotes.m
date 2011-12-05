% I = newton_cotes(f, a,b, degree, num_points)
%
% Integrates the function f on the interval [a,b] using the Newton-Cotes
% formula of the specified degree and (approximately) num_points evaluation
% points.
%
function I = newton_cotes(f, a,b, d, n)
n = ceil((n-d-1)/d);
w = d * (vander((0:d)/d)' \ (1./fliplr(1:d+1))')';
w = [w(1:end-1) repmat(w(1:end-1) + [w(end) zeros(1, d-1)], 1, n) w(end)];
x = linspace(a,b, length(w));
h = x(2) - x(1);
I = h * w * f(x)';
