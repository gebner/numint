function [I,Err,A]=mulDimNumInt(f,a,b,epsilon,w,p,alpha)
if nargin<7
    [p,alpha] = gauss_arrays(2, 12);
    w = newton_cotes_weights(2, 12);
end
s=size(a);
if s==[1,1]
    [I,Err,A]=NumIntMul(@(x) dimOne(f,x),a,b,epsilon,w,p,alpha,4,6,4,8);
else
    [I,Err,A]=NumIntMul(@(x) dimDown(f,x,a(2:end),b(2:end),epsilon,w,p,alpha),a(1),b(1),epsilon,w,p,alpha,4,6,4,8);
end

function [z,A]=dimDown(f,X,a,b,epsilon,w,p,alpha)
z=[X;X];
n=1;
A=zeros(2,4);
for x=X
    [z(1,n),z(2,n),At]=mulDimNumInt(@(y) f([x,y]),a,b,epsilon,w,p,alpha);
    A=A+At;
    n=n+1;
end

function [z,A]=dimOne(f,X)
z=X;
A=zeros(2,4);
n=1;
for x=X
    z(n)=f(x);
    n=n+1;
end