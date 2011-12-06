function [I, failure] = newton_cotes(f, a, b, epsilon, w)
% NEWTON_COTES evaluates the integral \int_a^b f(x) dx using the newton cotes method
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

iMin = 10;
iMax = 12;

finished = 0;
val = 0;
i = iMin;

while (~finished) % degree for fitting
   oldval = val;
   val = newton_cotes_eval(f, a, b, i, w);
    
    %get failure estimation
   failure  = 2*abs((val -oldval) / val);  
   if (i < iMin)
       finished = 0;
   elseif(failure < epsilon)
       finished = 1;
   elseif (i > iMax)
       finished = 1;
   end
   i = i+1;
end

I = val;
end

function I = newton_cotes_eval(f, a,b, d, w)
% NEWTON_COTES_EVAL Integrates the function f on the interval [a,b] using the Newton-Cotes
% formula of the specified degree on n intervals
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] d          degree for interpolation
% @param[in] n          number of intervals
% @param[out] I         integral value

% evaluate weights
%w = d * (vander((0:d)/d)' \ (1./fliplr(1:d+1))')'; %solve LGS
%w = [w(1:end-1) (w(1:end-1) + [w(end) zeros(1, d-1)]) w(end)];
w1 = w(d, (end-2*d):end);
% nodes
x = linspace(a,b, length(w1));
h = x(2) - x(1);% step width

%evaluate integral 
I = h * w1 * f(x)';

end
