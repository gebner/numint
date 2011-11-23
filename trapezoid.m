function I = trapezoid(f, a,b, n),
I = newton_cotes(f, a,b, trapezoid_weights(n));

function weights = trapezoid_weights(n),
weights = [1/2 ones(1, n-2) 1/2];
