function [I,Err,A]=NumInt(f,a,b,epsilon)
% NUMINT performs numeric integration
% @param[in]    f       function pointer
% @param[in]    a       lower bound
% @param[in]    b       upper bound
% @param[in]    epsilon accuracy
% @param[out]   I       Integral value
% @param[out]   Err     error (a posteriori)
% @param[out] A     statistics array

    tic;    
    
    Gmin = 8; Gmax = 10; %index constants for gauss integral
    Rmin = 6; Rmax = 8; %index constants for romberg integral
    Nmin = 8; Nmax = 12; %index constants for newton cotes integral
    
    if (epsilon<10^-15)
        epsilon=10^-15;
    end
    
    %eval weigths
    [x,alpha] = gauss_arrays(Gmin, Gmax);
    w = newton_cotes_weights(Nmin, Nmax);
    
    %perform evaluation
    [I,Err,A]=NumIntStep(f,a,b,epsilon,x,alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
    I=sum(I); % sum of all parts
    
    toc;
end

function [I,Err,A]=NumIntStep(f,a,b,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax)
%  NUMINTSTEP divides the intervall and chooses the algorithm
% @param[in]    f       function pointer
% @param[in]    a       lower bound
% @param[in]    b       upper bound
% @param[in]    epsilon accuracy
% @param[in]    x       matrix of legendre roots
% @param[in]    alpha   matrix of legendre weights
% @param[in]    w       weights for newton cotes
% @param[in]    Nmin       minimum index for newton cotes
% @param[in]    Nmax       maximum index for newton cotes
% @param[in]    Rmin       minimum index for romberg
% @param[in]    Rmax       maximum index for romberg
% @param[in]    Gmin       minimum index for gauss
% @param[in]    Gmax       maximum index for gauss
% @param[out]    I      array containing the integrals of all intervall parts
% @param[out]    Err    error (a posteriori)
% @param[out]   A     statistics array

 [R,eR]=romberg(f,a,b,epsilon,Rmin, Rmax); %try romberg
 if (isnan(eR) || isinf(eR))
     eR=inf;
 end

 [N, eN] = newton_cotes(f, a, b, epsilon, w, Nmin, Nmax); %try newton cotes
 if isnan(eN) || isinf(eN)
     eN=inf;
 end
 
 [G,eG]=gauss(f,a,b,epsilon,x,alpha, Gmin, Gmax); %try gauss
 if isnan(eG) || isinf(eG)
     eG=inf;
 end
 
 % decide what to do
 if (min([eR,eN,eG])<epsilon || abs(b-a) < epsilon) % is epsilon plausible
     if (~isnan(R) && ~isinf(R) && eR<=eN && eR<=eG) % romberg integration worked best
         if (((isnan(N)||abs(R-N)<=abs(N*eN)) && (isnan(G)||abs(R-G)<=abs(G*eG))) || abs(b-a) < epsilon) % romberg result is in error ranges
             I=R;
             Err=eR;
             A=zeros(2,4); % statistics array
             A(2,1)=A(2,1)+1;
         else % romberg not in error ranges -> somethings wrong -> divide
             [I1,E1,A1]=NumIntStep(f,a,(a+b)/2,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             [I2,E2,A2]=NumIntStep(f,(a+b)/2,b,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             I=[I1,I2];
             A = A1 + A2;
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             A(1,1)=A(1,1)+1;
         end
     elseif (~isnan(N) && ~isinf(N) && eN<=eG && eN<=eR) % newton cotes worked best
         if (((isnan(G)||abs(N-G)<=abs(G*eG)) && (isnan(R)||abs(R-N)<=abs(R*eR))) || abs(b-a) < epsilon) % newton cotes result is in error ranges
            I=N;
            Err=eN;
            A=zeros(2,4); % statistics array
            A(2,2)=A(2,2)+1;
         else % newton cotes not in error ranges -> somethings wrong -> divide
             [I1,E1,A1]=NumIntStep(f,a,(a+b)/2,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             [I2,E2,A2]=NumIntStep(f,(a+b)/2,b,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             I=[I1,I2];
             A = A1 + A2;
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             A(1,2)=A(1,2)+1;
         end
     elseif (~isnan(G) && ~isinf(G)) % gauss worked best
        if (((isnan(R)||abs(R-G)<=abs(R*eR)) && (isnan(N)||abs(N-G)<=abs(N*eN))) || abs(b-a)< epsilon) % gauss result is in error ranges
            I=G;
            Err=eG;
            A=zeros(2,4); % statistics array
            A(2,3)=A(2,3)+1;
         else % gauss not in error ranges -> somethings wrong -> divide
             [I1,E1,A1]=NumIntStep(f,a,(a+b)/2,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             [I2,E2,A2]=NumIntStep(f,(a+b)/2,b,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
             I=[I1,I2];
             A = A1 + A2;
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             A(1,3)=A(1,3)+1;
         end
     end
 else % epsilon too small
     [I1,E1,A1]=NumIntStep(f,a,(a+b)/2,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
     [I2,E2,A2]=NumIntStep(f,(a+b)/2,b,epsilon, x, alpha, w, Nmin, Nmax, Rmin, Rmax, Gmin, Gmax);
     I=[I1,I2];
     A = A1 + A2;
     A(1,4)=A(1,4)+1;
     int1 = sum(I1);
     int2 = sum(I2);
     Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
 end
end 
