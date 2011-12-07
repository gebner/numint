function [integral , failure] = romberg(f, a, b, epsilon, iMin, iMax)
%ROMBERG evaluates the integral \int_a^b f(x) dx using the romberg method
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[in] iMin       minimum index (optional)
% @param[in] iMax       maximum index (optional)
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

if(nargin <= 6)
    iMin = 1; % min index
    iMax = 8; % max index
end

failure = Inf;

d = zeros(iMax, 1); % vector for diagonal values
d(1) = Inf;
i = 0;
while ((i < iMin) ||(failure > epsilon && i <= iMax))
   
   % evaluate area by summing trapezes
   d(i+1) = evalTrapezArea(f, a, b, 2^i);
   
   prev = d(1);
   % refresh previous last elements vector
   for j = i:(-1):1
     d(j) = d(j+1) + (d(j+1)-d(j))/(2^(2*(i-j+1))-1);
   end
   
   % get failure estimation
   failure = abs((d(1) - prev)/prev);
   
   i = i+1; % new stepwidth = (b-a)/2^i
end

integral = d(1);
     
end


function area = evalTrapezArea(f, a, b, i)
% EVALTRAPEZAREA evaluates the area of sums of i trapezes
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] i          number of trapezes
% @param[out]  area     area of trapezes

   % evaluate i-th trapezsum
   trapezsum = f(a) + f(b); %sum of trapez heights
   for j = 1:(i-1)
       pos = a + j*(b-a)/i; % position
       trapezsum = trapezsum + 2*f(pos);
   end
   area = (b-a)/(2*i) * trapezsum; 
end

