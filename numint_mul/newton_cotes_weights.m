function [ w ] = newton_cotes_weights(dMin, dMax)
%NEWTON_COTES_WEIGHTS derives the weights for the newton-cotes integration
%and stores them in a matrix
% weights for degree i are stored in the i-th row in the last 2*i+1
% elements
% @param[in]    dMin        minimal degree
% @param[in]    dMax        maximal degree
% @param[out]   w           matrix containing the weights, where the
w = zeros(dMax, 2*dMax+1);

for i = dMin:dMax
    w_help = i * (vander((0:i)/i)' \ (1./fliplr(1:i+1))')'; %solve LGS
    w(i,(2*dMax+1-2*i):2*dMax+1) = [w_help(1:end-1) (w_help(1:end-1) + [w_help(end) zeros(1, i-1)]) w_help(end)];
end

end

