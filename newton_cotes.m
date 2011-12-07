function [val, failure] = newton_cotes(f, a, b, epsilon, w, iMin, iMax)
% NEWTON_COTES evaluates the integral \int_a^b f(x) dx using the newton cotes method
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[in] w          array with coefficients (optional)
% @param[in] iMin       minimum index (optional)
% @param[in] iMax       maximum index (optional)
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

if (nargin < 7)
    iMin = 8;
    iMax = 12;
end
if (nargin == 4)
    w = newton_cotes_weight(iMin, iMax);
end

failure = Inf;
val = Inf;

i = iMin;
while (failure > epsilon && i <= iMax) % degree for fitting
   oldval = val;
   %eval integral
   val = newton_cotes_eval(f, a, b, i, w);
    
    %get failure estimation
   failure  = abs((val - oldval) / val);
   
   i = i+2;
end

end

function I = newton_cotes_eval(f, a, b, d, w)
% NEWTON_COTES_EVAL Integrates the function f on the interval [a,b] using the Newton-Cotes
% formula of the specified degree
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] d          degree for interpolation
% @param[in] w          array of coefficents
% @param[out] I         integral value

% evaluate weights
w1 = w(d, (end-2*d):end);

% nodes
x = linspace(a,b, length(w1));
h = x(2) - x(1);% step width

%evaluate integral 
I = h * w1 * f(x)';

end
