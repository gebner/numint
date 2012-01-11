
I=[];
E=[];
F=[];
T=[];
x=2:15;
y=4.^(-x);
for eps=y
   tic;
   [i,e]=NumInt(@(x) sin(1./x),0,1,eps);
   t=toc;
   f=abs((0.5040670619069283719898561177411482-i)/i);
   I=[I,i];
   E=[E,e];
   F=[F,f];
   T=[T,t];
end
