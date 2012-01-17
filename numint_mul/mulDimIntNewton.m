function [I,Err,A]=mulDimIntNewton(f,a,b,epsilon,A, Nmin, Nmax, w)
% MULDIMINTOPT multidimensional integration with preevaluated coefficients
% optimized version of mulDimInt
% @param[in]  f       function pointer
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy
% @param[in] Nmin       newton cotes minimum index (optional, if preevaluated)
% @param[in] Nmax       newton cotes maximum index (optional, if preevaluated)
% @param[in] w          newton cotes weights, (optional) preevaluated for
% optimization
% @param[out] I       value of integral
% @param[out] Err     error approximation
% @param[in,out] A       statistics array (optional)

if (nargin <5)
    A = zeros(2,4);
end
if (nargin<8)
    Nmin = 4; Nmax = 8;
    w = newton_cotes_weights(Nmin, Nmax);
end
s=length(a);
if (s ~= length(b)) %number of lower != number of upper bounds
   error('Boundries dont match!');
end
if s==1 %one dimensional
    [I,Err,A]=newton_cotes(@(x) mapWithStats(f,x),a,b,epsilon,A,Nmin,Nmax, w);
else %more than one dimension left
    [I,Err,A]=newton_cotes(@(x) dimDown(f,x,a(2:end),b(2:end),epsilon, A,Nmin, Nmax,w),a(1),b(1),epsilon,A, Nmin,Nmax,w);
end 
end

% DIMDOWN fixes the first variable and integrates by fubini
% @param[in]  f       function pointer
% @param[in]  X       array
% @param[in] Int      function pointer to integrator (working also on arrays)
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy
% @param[in] Nmin       newton cotes minimum index 
% @param[in] Nmax       newton cotes maximum index 
% @param[in] w          newton cotes weights, preevaluated for optimization
% @param[in,out] A    statistics array 
% @param[out] z       2 x length(X) array containing integral values and errors
function [z,A] =dimDown(f,X,a,b,epsilon,A,Nmin, Nmax, w)

len = size(X,2);
z=zeros(2,len);

for n = 1:len % evaluate for fixed first variable x = x(n)
    [z(1,n), z(2,n), A] = mulDimIntNewton(@(y) f([X(n),y]),a,b,epsilon, A,Nmin, Nmax, w);
end

end
