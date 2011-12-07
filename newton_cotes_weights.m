function [ w ] = newton_cotes_weights(dMin, dMax )
%NEWTON_COTES_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
w = zeros(dMax, 2*dMax+1);

for i = dMin:dMax
    w_help = i * (vander((0:i)/i)' \ (1./fliplr(1:i+1))')'; %solve LGS
    w(i,(2*dMax+1-2*i):2*dMax+1) = [w_help(1:end-1) (w_help(1:end-1) + [w_help(end) zeros(1, i-1)]) w_help(end)];
end

end

