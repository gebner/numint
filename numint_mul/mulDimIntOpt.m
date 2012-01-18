function [I,Err,A]=mulDimIntOpt(f,a,b,epsilon,A,Gmin, Gmax, Nmin, Nmax, w,p,alpha)
% MULDIMINTOPT multidimensional integration with preevaluated coefficients
% optimized version of mulDimInt
% @param[in]  f       function pointer
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy
% @param[in] Gmin       gauss minimum index (optional, if preevaluated)
% @param[in] Gmax       gauss maximum index (optional, if preevaluated)
% @param[in] Nmin       newton cotes minimum index (optional, if preevaluated)
% @param[in] Nmax       newton cotes maximum index (optional, if preevaluated)
% @param[in] w          newton cotes weights, (optional) preevaluated for optimization
% @param[in] p          gauss roots, (optional) preevaluated for
% optimization
% @param[in] alpha      gauss weights, (optional) preevaluated for optimization
% @param[out] I       value of integral
% @param[out] Err     error approximation
% @param[in,out] A       statistics array (optional)

if (nargin <5)
    A = zeros(2,4);
end
if (nargin<12)
    Gmin = 2;Gmax = 4;
    [p,alpha] = gauss_arrays(Gmin, Gmax+1);
    Nmin = 2; Nmax = 6;
    w = newton_cotes_weights(Nmin, Nmax);
end
s=length(a);
if (s ~= length(b)) %number of lower != number of upper bounds
   error('Boundries dont match!');
end
if s==1 %one dimensional
    [I,Err,A]=NumInt(@(x) mapWithStats(f,x),a,b,epsilon,A,Gmin,Gmax,Nmin,Nmax, w,p,alpha);
else %more than one dimension left
    [I,Err,A]=NumInt(@(x) dimDown(f,x,a(2:end),b(2:end),epsilon, A, Gmin, Gmax, Nmin, Nmax,w,p,alpha),a(1),b(1),epsilon,A, Gmin,Gmax,Nmin,Nmax,w,p,alpha);
end 
end

% DIMDOWN fixes the first variable and integrates by fubini
% @param[in]  f       function pointer
% @param[in]  X       array
% @param[in] Int      function pointer to integrator (working also on arrays)
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy
% @param[in] Gmin       gauss minimum index 
% @param[in] Gmax       gauss maximum index 
% @param[in] Nmin       newton cotes minimum index 
% @param[in] Nmax       newton cotes maximum index 
% @param[in] w          newton cotes weights, preevaluated for optimization
% @param[in] p          gauss roots, preevaluated for optimization
% @param[in] alpha      gauss weights, preevaluated for optimization
% @param[in,out] A    statistics array 
% @param[out] z       2 x length(X) array containing integral values and errors
function [z,A] =dimDown(f,X,a,b,epsilon,A,Gmin, Gmax, Nmin, Nmax, w, p, alpha)

len = size(X,2);
z=zeros(2,len);

for n = 1:len % evaluate for fixed first variable x = x(n)
    [z(1,n), z(2,n), A] = mulDimIntOpt(@(y) f([X(n),y]),a,b,epsilon, A,Gmin, Gmax, Nmin, Nmax, w, p, alpha);
end

end
