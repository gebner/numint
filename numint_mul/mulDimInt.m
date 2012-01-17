function [I,Err, A]=mulDimInt(f,Integrator,a,b,epsilon, A)
% MULDIMINT perform multidimensional integration by using fubini
% @param[in] f        function pointer
% @param[in] Int      function pointer to integrator (working also on arrays)
%                     has to be of form @(f,a,b,eps,A) Integrator(f,a,b,eps,A)
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy   
% @param[in,out] A    statistics array  
% @param[out] I       value of integral (optional)
% @param[out] Err     error approximation
if (nargin < 6)
    A = zeros(2,4);
end

s=length(a);
if (s ~= length(b)) %number of lower != number of upper bounds
   error('Boundries dont match!');
end
   
if (s==1) % one dimensional integral
    [I,Err,A]=Integrator(@(x) mapWithStats(f,x),a,b,epsilon, A); %arrayfun maps f on array x
else % more than one dimensional integral
    [I,Err,A]=Integrator(@(x) dimDown(f,x,Integrator,a(2:end),b(2:end),epsilon, A),a(1),b(1),epsilon, A);
end

end

function [z, A] = dimDown(f,X,Integrator,a,b,epsilon, A)
% DIMDOWN fixes the first variable and integrates by fubini
% @param[in]  f       function pointer
% @param[in]  X       array
% @param[in] Int      function pointer to integrator (working also on arrays)
% @param[in] a        lower bound array
% @param[in] b        upper bound array
% @param[in] epsilon  accuracy   
% @param[in,out] A    statistics array  
% @param[out] z       2 x length(X) array containing integral values and
% errors

len = size(X,2);
z = zeros(2, len);

for n = 1:len % evaluate for fixed first variable x = X(n)
    [z(1,n), z(2,n),A] = mulDimInt(@(y) f([X(n),y]),Integrator,a,b,epsilon, A);
end

end


