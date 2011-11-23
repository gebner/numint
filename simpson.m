function I = simpson(f, a,b, n),
I = newton_cotes(f, a,b, simpson_weights(n));

function weights = simpson_weights(n),
weights = [1/3 4/3 repmat([2/3 4/3], 1, ceil((n-3)/2)) 1/3];
