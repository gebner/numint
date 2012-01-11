function [I,Err]=mulDimInt(f,Int,a,b,epsilon)
s=size(a);
if s==[1,1]
    [I,Err]=Int(@(x) dimOne(f,x),a,b,epsilon);
else
    [I,Err]=Int(@(x) dimDown(f,x,Int,a(2:end),b(2:end),epsilon),a(1),b(1),epsilon);
end



function z=dimOne(f,X)
z=X;
n=1;
for x=X
    z(n)=f(x);
    n=n+1;
end

function z=dimDown(f,X,Int,a,b,epsilon)
z=[X;X];
n=1;
for x=X
    [z(1,n),z(2,n)]=mulDimInt(@(y) f([x,y]),Int,a,b,epsilon);
    n=n+1;
end

