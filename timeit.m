% [avg, n] = timeit(handle)
%
% Runs handle several times and returns the average run time in seconds.
%
% The number of times the handle is run is chosen to exceed 50ms in total and
% is returned as n.
% 
% > [avg, n] = timeit(@() vander(1:10) \ (1:10)')
% avg =  0.0002247
% n =  1000
%
function [avg, n] = timeit(handle),
for n = [3 1e1 1e2 1e3 1e4 1e5],
  start = time;
  for i = 1:n,
    handle();
  end
  stop = time;

  dur = stop - start;
  avg = dur / n;

  if dur >= 0.05,
    break
  end
end
