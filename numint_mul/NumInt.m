function [I,Err,A]=NumInt(f,a,b,epsilon,A,Gmin,Gmax,Nmin,Nmax,w,x,alpha)
% NUMINT performs numeric integration
% @param[in]    f       function pointer
% @param[in]    a       lower bound
% @param[in]    b       upper bound
% @param[in]    epsilon accuracy
% @param[out]   I       Integral value
% @param[out]   Err     error (a posteriori)
% @param[in,out] A     statistics array (optional)
% @param[in] w          newton cotes weights, (optional) preevaluated for optimization
% @param[in] x          gauss roots, (optional) preevaluated for optimization
% @param[in] alpha      gauss weights, (optional) preevaluated for optimization
% @param[in] Gmin       gauss minimum index (optional, if preevaluated)
% @param[in] Gmax       gauss maximum index (optional, if preevaluated)
% @param[in] Nmin       newton cotes minimum index (optional, if preevaluated)
% @param[in] Nmax       newton cotes maximum index (optional, if preevaluated)
  
    if(nargin<5)
      A = zeros(2,4);
      f = @(x) mapWithStats(f,x); % add dummy stats array
    end
    %set boundries
    if nargin<7
      Gmin = 6; Gmax = 8; %index constants for gauss integral
    end
    if nargin<9
      Nmin = 6; Nmax = 10; %index constants for newton cotes integral
    end
    %eval weigths
    if nargin<12
        [x,alpha] = gauss_arrays(Gmin, Gmax+1);
        w = newton_cotes_weights(Nmin, Nmax);
    end
    if (epsilon<10^-15)
        epsilon=10^-15;
    end
    
    %perform evaluation
    [I,Err,A]=NumIntStep(f,a,b,epsilon,A,x,alpha, w, Nmin, Nmax, Gmin, Gmax);
    I=sum(I); % sum of all parts
    
end

function [I,Err,A]=NumIntStep(f,a,b,epsilon, A, x, alpha, w, Nmin, Nmax, Gmin, Gmax)
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

 [N, eN, A] = newton_cotes(f, a, b, epsilon,A, Nmin, Nmax, w); %try newton cotes
 if isnan(eN) || isinf(eN)
     eN=inf;
 end
 [G,eG,A]=gauss(f,a,b,epsilon,A, Gmin, Gmax, x, alpha); %try gauss
 if isnan(eG) || isinf(eG)
     eG=inf;
 end
 % decide what to do
 if ((min([eN,eG])<epsilon && ((~isnan(G) && ~isinf(G))||(~isnan(N) && ~isinf(N))))|| abs(b-a) < epsilon) % is epsilon plausible
     if (~isnan(N) && ~isinf(N) && (eN<=eG || (isnan(eG) || isinf(eG)))) % newton cotes worked best
         if ((((isnan(G)||abs(N-G)<=abs(G*eG)))) || abs(b-a) < epsilon) % newton cotes result is in error ranges
            I=N;
            Err=eN;
            A(2,2)=A(2,2)+1; % statistics array
         else % newton cotes not in error ranges -> somethings wrong -> divide
             [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
             [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
             I=[I1,I2];
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             A(1,2)=A(1,2)+1;
         end
     elseif (~isnan(G) && ~isinf(G)) % gauss worked best
        if (((isnan(N)||abs(N-G)<=abs(N*eN))) || abs(b-a)< epsilon) % gauss result is in error ranges
            I=G;
            Err=eG;
            A(2,3)=A(2,3)+1;% statistics array
         else % gauss not in error ranges -> somethings wrong -> divide
             [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
             [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
             I=[I1,I2];
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             A(1,3)=A(1,3)+1;
        end
     else
         error('no valid result using the implemented integration methods');
     end
 else % epsilon too small
     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w, max(Nmax-2,Nmin), Nmax, max(Gmax-2,Gmin), Gmax);
     I=[I1,I2];
     A(1,4)=A(1,4)+1;
     int1 = sum(I1);
     int2 = sum(I2);
     Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
 end
end 