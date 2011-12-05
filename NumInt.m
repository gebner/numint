function [I,Err]=NumInt(f,a,b,epsilon)
    t=tic;
    if epsilon<10^-15
        epsilon=10^-15;
    end
    A=zeros(2,4);
    [I,Err,A]=NumIntStep(f,a,b,epsilon,A);
    size(I)
    I=sum(I);
    tic-t
    A
end

function [I,Err,A]=NumIntStep(f,a,b,epsilon,A)
%     if abs(b-a)<max(epsilon,10^-15)
%          if abs(f(a))<=abs(f(b))
%              I(i)=(b-a)*f(a);
%          else
%              I(i)=(b-a)*f(b);
%          end
%     Err=epsilon;
%     i=i+1;
%     A(2,4)=A(2,4)+1;
%     else
         [R,eR]=romberg(f,a,b,epsilon);
         if isnan(eR) || isinf(eR)
             eR=inf;
         end
         N6=newton_cotes(f,a,b,10,7);
%         N7=newton_cotes(f,a,b,11,8);
         N=newton_cotes(f,a,b,12,9);
         eN=abs(abs(N-N6)/N);
         if isnan(eN) || isinf(eN)
             eN=inf;
         end
         [G,eG]=integrieren_gauss(f,a,b,epsilon);
         if isnan(eG) || isinf(eG)
             eG=inf;
         end
         if min([eR,eN,eG])<epsilon || abs(b-a)<max(epsilon,10^-15)
             if ~isnan(R) && ~isinf(R) && eR<=eN && eR<=eG
                 if (abs(R-N)<=abs(N*eN) && abs(R-G)<=abs(G*eG)) || abs(b-a)<max(epsilon,10^-15)
                    I=R;
                    Err=eR;
                    A(2,1)=A(2,1)+1;
                 else
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A);
                     I=[I1,I2];
                     Err=max(E1,E2);
                     A(1,1)=A(1,1)+1;
%                     [a,b,1;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             else
             if ~isnan(N) && ~isinf(N) && eN<=eG && eN<=eR
                 if (abs(N-G)<=abs(G*eG) && abs(R-N)<=abs(R*eR)) || abs(b-a)<max(epsilon,10^-15)
                    I=N;
                    Err=eN;
                    A(2,2)=A(2,2)+1;
                 else
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A);
                     I=[I1,I2];
                     Err=max(E1,E2);
                     A(1,2)=A(1,2)+1;
%                     [a,b,2;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             else
                 if (abs(R-G)<=abs(R*eR) && abs(N-G)<=abs(N*eN)) || abs(b-a)<max(epsilon,10^-15)
                    I=G;
                    Err=eG;
                    A(2,3)=A(2,3)+1;
                 else
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A);
                     I=[I1,I2];
                     Err=max(E1,E2);
                     A(1,3)=A(1,3)+1;
%                     [a,b,3;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             end
             end
         else
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A);
                     I=[I1,I2];
                     A(1,4)=A(1,4)+1;
                     Err=max(E1,E2);
%                     [a,b,4;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
         end
%     end
end 