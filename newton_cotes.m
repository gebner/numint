function I = newton_cotes(f, a,b, weights),
x = linspace(a,b, length(weights));
h = x(2) - x(1);
I = h * weights * f(x)';
