function [val, failure,A] = newton_cotes(f, a, b, epsilon, iMin, iMax, w)
% NEWTON_COTES evaluates the integral \int_a^b f(x) dx using the newton cotes method
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[in] w          array with coefficients (optional)
% @param[in] iMin       minimum index (optional)
% @param[in] iMax       maximum index (optional)
% @param[in,out] A     statistics array (optional)
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

A=zeros(2,4);
if (nargin < 5)
    f = @(x) mapWithStats(f,x); % add dummy stats array
end

if (nargin < 6)
    iMin = 8;
    iMax = 12;
end
if (nargin < 7)
    w = newton_cotes_weights(iMin, iMax);
end

failure = Inf;
val = Inf;
i = iMin;
while (failure > epsilon && i <= iMax) % degree for fitting
   oldval = val;
   %eval integral
   [val,Err,At] = newton_cotes_eval(f, a, b, i, w);
   A = A+At; %refresh statistics
   
    %get failure estimation
    if(val~=0)
        failure  = abs((abs(val - oldval)+Err) / val);
    else
        failure = abs(abs(val-oldval)+Err);
    end
   
   i = i+2;
end

end

function [I,Err,A] = newton_cotes_eval(f, a, b, d, w)
% NEWTON_COTES_EVAL Integrates the function f on the interval [a,b] using the Newton-Cotes
% formula of the specified degree
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] d          degree for interpolation
% @param[in] w          array of coefficents
% @param[out] I         integral value
% @param[out] Err       error approximation
% @param[out] A         statistics array

% evaluate weights
w1 = w(d, (end-2*d):end);

% nodes
x = linspace(a,b, length(w1));
h = x(2) - x(1);% step width
%evaluate integral 
[I,A]=f(x);
s=size(I);

if s(1)==2
    Err=abs(I(2,:).*I(1,:));
    Err=h*abs(w1)*Err';
else
    Err=0;
end
I = h * w1 * I(1,:)';
end
