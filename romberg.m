function [integral , failure] = romberg(f, a, b, epsilon)
%ROMBERG evaluates the integral \int_a^b f(x) dx using the romberg method
% @param[in] f          function pointer
% @param[in] a          lower bound
% @param[in] b          upper bound
% @param[in] epsilon    accuracy
% @param[out] integral  evaluated value
% @param[out] failure   evaluation failure (a posteriori)

delta = 1/10^9; % delta value for undefined positions
iMin = 8; % min index
iMax = 10; % max index

finished = 0; %indicator for end of recursion
i = 1;

d = zeros(iMin, 1); %diagonal values
d(1) = (b-a)/2 * (f(a) + f(b));
while (finished ~= 1)
   i = i+1; % new stepwidth = (b-a)/2^(i-1)
   
   % evaluate i-th trapezsum
   trapezsum = f(a) + f(b); %sum of trapez heights
   pos = a;
   for j = 1:(2^(i-1)-1)
       pos = pos + (b-a)/2^(i-1); % position
       trapezsum = trapezsum + 2*f(pos);
   end
   d(i) = (b-a)/2^i * trapezsum; 
   
   prev = d(1);
   % refresh previous last elements vector
   for j = (i-1):(-1):1
     d(j) = d(j+1) + (d(j+1)-d(j))/(2^(2*(i-j))-1);
   end
   
   % get failure estimation
   failure = abs((d(1)-prev)/prev);
   if (i < iMin)
       finished = 0;
   elseif(failure < epsilon)
       finished = 1;
   elseif (i > iMax)
       finished = 1;
   end
end

integral = d(1);
     
end

