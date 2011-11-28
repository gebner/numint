function [integral , failure] = romberg(f, a, b, epsilon)
%ROMBERG Summary of this function goes here
%   Detailed explanation goes here


delta = 1/10^9;
iMin = 3; % min index
d = zeros(iMin, 1); %diagonal values
iMax = 22; % max index
finished = 0;
i = 1;

if (isnan(f(a)) || isinf(f(a)))
   a = a + (b-a)*delta;
end
if (isnan(f(b)) || isinf(f(b)))
   b = b - (b-a)*delta;
end   

d(1) = (b-a)/2 * (f(a) + f(b));
while (finished ~= 1)
   i = i+1; % new stepwidth = (b-a)/2^(i-1)
   
   % evaluate i-th trapezsum
   trapezsum = f(a) + f(b); %sum of trapez heights
   pos = a;
   for j = 1:(2^(i-1)-1)
       pos = pos + (b-a)/2^(i-1); % position
       fpos = 2*f(pos);
       if (isnan(fpos) || isinf(fpos)) 
           fpos = f(pos + 2^(-i-1)) + f(pos - 2^(-i-1));
       end
       trapezsum = trapezsum + fpos;
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

