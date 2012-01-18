function [X, A] = mapWithStats(f,X)
%MAPWITHSTATS map with dummy statistic array
X = arrayfun(f,X);
A = zeros(2,4);

end

