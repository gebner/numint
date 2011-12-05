function [I,Err]=NumInt2(f,a,b,epsilon)
    t=tic;
    if epsilon<10^-15
        epsilon=10^-15;
    end
    [I,Err,tN,tR,tG]=NumIntStep(f,a,b,epsilon);
    tic-t
    tN=tN
    tR=tR
    tG=tG
end

function [I,Err,tN,tR,tG]=NumIntStep(f,a,b,epsilon)
     if abs(b-a)<max(epsilon,10^-15)
          if abs(f(a))<=abs(f(b))
              I=(b-a)*f(a);
          else
              I=(b-a)*f(b);
          end
        Err=epsilon;
        tN=0;
        tG=0;
        tR=0;
     else
         tic;
         [R,eR]=romberg(f,a,b,epsilon);
         tR=toc;
         tic;
         N6=newton_cotes(f,a,b,6,7);
         N7=newton_cotes(f,a,b,7,8);
         N=newton_cotes(f,a,b,8,9);
         eN=abs(min(abs(N-N7),abs(N-N6))/N);
         tN=toc;
         tic;
         [G,eG]=integrieren_gauss(f,a,b,epsilon);
         tG=toc;
         if min([eR,eN,eG])<epsilon
             if eR<=eN && eR<=eG
                 if abs(R-N)<=eN && abs(R-G)<=eN
                    I=R;
                    Err=eR;
                 else
                     [I1,E1,tN1,tR1,tG1]=NumIntStep(f,a,(a+b)/2,epsilon);
                     [I2,E2,tN2,tR2,tG2]=NumIntStep(f,(a+b)/2,b,epsilon);
                     I=I1+I2;
                     Err=max(E1,E2);
                     tN=tN+tN1+tN2;
                     tR=tR+tR1+tR2;
                     tG=tG+tG1+tG2;
                 end
             else
             if eR<=eG && eR<=eN
                 if abs(R-G)<=eR && abs(R-N)<=eR
                    I=N;
                    Err=eN;
                 else
                     
                     [I1,E1,tN1,tR1,tG1]=NumIntStep(f,a,(a+b)/2,epsilon);
                     [I2,E2,tN2,tR2,tG2]=NumIntStep(f,(a+b)/2,b,epsilon);
                     I=I1+I2;
                     Err=max(E1,E2);
                     tN=tN+tN1+tN2;
                     tR=tR+tR1+tR2;
                     tG=tG+tG1+tG2;
                 end
             else
                 if abs(R-G)<=eG && abs(N-G)<=eG
                    I=G;
                    Err=eG;
                 else
                     [I1,E1,tN1,tR1,tG1]=NumIntStep(f,a,(a+b)/2,epsilon);
                     [I2,E2,tN2,tR2,tG2]=NumIntStep(f,(a+b)/2,b,epsilon);
                     I=I1+I2;
                     Err=max(E1,E2);
                     tN=tN+tN1+tN2;
                     tR=tR+tR1+tR2;
                     tG=tG+tG1+tG2;
                 end
             end
             end
         else
                     [I1,E1,tN1,tR1,tG1]=NumIntStep(f,a,(a+b)/2,epsilon);
                     [I2,E2,tN2,tR2,tG2]=NumIntStep(f,(a+b)/2,b,epsilon);
                     I=I1+I2;
                     Err=max(E1,E2);
                     tN=tN+tN1+tN2;
                     tR=tR+tR1+tR2;
                     tG=tG+tG1+tG2;
         end
     end
end 
             