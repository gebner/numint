function [I,Err,A]=NumInt(f,a,b,epsilon)
% NUMINT performs numeric integration
% @param[in]    f       function pointer
% @param[in]    a       lower bound
% @param[in]    b       upper bound
% @param[in]    epsilon accuracy
% @param[out]   I       Integral value
% @param[out]   Err     error (a posteriori)
% @param[in, out] A     statistics array
    tic;
    if (epsilon<10^-15)
        epsilon=10^-15;
    end
    [x,alpha] = legendre_all(13);
    w = newton_cotes_weights(13);
    A=zeros(2,4); % statistics array
    [I,Err,A]=NumIntStep(f,a,b,epsilon,A,x,alpha, w);
    %size(I) % array containing the parts
    I=sum(I); % sum of all parts
    toc;
end

function [I,Err,A]=NumIntStep(f,a,b,epsilon,A, x, alpha, w)
%    NUMINTSTEP divides the intervall and chooses the algorithm
% @param[in]    f       function pointer
% @param[in]    a       lower bound
% @param[in]    b       upper bound
% @param[in]    epsilon accuracy
% @param[in, out] A     statistics array
% @param[in]    x       matrix of legendre roots
% @param[in]    alpha   matrix of legendre weights
% @param[out]    I      array containing the integrals of all intervall parts
% @param[out]    Err    error (a posteriori)

% old version
%     if (abs(b-a)< epsilon)
%          if abs(f(a))<=abs(f(b))
%              I(i)=(b-a)*f(a);
%          else
%              I(i)=(b-a)*f(b);
%          end
%     Err=epsilon;
%     i=i+1;
%     A(2,4)=A(2,4)+1;
%     else
         [R,eR]=romberg(f,a,b,epsilon); %try romberg
         if (isnan(eR) || isinf(eR))
             eR=inf;
         end
         
         [N, eN] = newton_cotes(f, a, b, epsilon, w); %try newton cotes
         if isnan(eN) || isinf(eN)
             eN=inf;
         end
         [G,eG]=gauss(f,a,b,epsilon,x,alpha); %try gauss
         if isnan(eG) || isinf(eG)
             eG=inf;
         end
         % decide what to do
         if (min([eR,eN,eG])<epsilon || abs(b-a) < epsilon) % is epsilon plausible
             if (~isnan(R) && ~isinf(R) && eR<=eN && eR<=eG) % romberg integration worked best
                 if ((  (isnan(N)||abs(R-N)<=abs(N*eN)) && (isnan(G)||abs(R-G)<=abs(G*eG))) || abs(b-a) < epsilon) % romberg result is in error ranges
                     I=R;
                     Err=eR;
                     A(2,1)=A(2,1)+1;
                 else % romberg not in error ranges -> somethings wrong -> divide
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w);
                     I=[I1,I2];
                     int1 = sum(I1);
                     int2 = sum(I2);
                     Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
                     %Err = max(E1, E2);
                     A(1,1)=A(1,1)+1;
%                     [a,b,1;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             elseif (~isnan(N) && ~isinf(N) && eN<=eG && eN<=eR) % newton cotes worked best
                 if (((isnan(G)||abs(N-G)<=abs(G*eG)) && (isnan(R)||abs(R-N)<=abs(R*eR))) || abs(b-a) < epsilon) % newton cotes result is in error ranges
                    I=N;
                    Err=eN;
                    A(2,2)=A(2,2)+1;
                 else % newton cotes not in error ranges -> somethings wrong -> divide
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w);
                     I=[I1,I2];
                     int1 = sum(I1);
                     int2 = sum(I2);
                     Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
                     %Err = max(E1, E2);
                     A(1,2)=A(1,2)+1;
%                     [a,b,2;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             elseif (~isnan(G) && ~isinf(G)) % gauss worked best
                if (((isnan(R)||abs(R-G)<=abs(R*eR)) && (isnan(N)||abs(N-G)<=abs(N*eN))) || abs(b-a)< epsilon) % gauss result is in error ranges
                    I=G;
                    Err=eG;
                    A(2,3)=A(2,3)+1;
                 else % gauss not in error ranges -> somethings wrong -> divide
                     [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w);
                     [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w);
                     I=[I1,I2];
                     int1 = sum(I1);
                     int2 = sum(I2);
                     Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
                     %Err = max(E1, E2);
                     A(1,3)=A(1,3)+1;
%                     [a,b,3;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
                 end
             end
         else % epsilon too small
             [I1,E1,A]=NumIntStep(f,a,(a+b)/2,epsilon,A, x, alpha, w);
             [I2,E2,A]=NumIntStep(f,(a+b)/2,b,epsilon,A, x, alpha, w);
             I=[I1,I2];
             A(1,4)=A(1,4)+1;
             int1 = sum(I1);
             int2 = sum(I2);
             Err= (abs(int1*E1) + abs(int2*E2))/abs(int1+int2);
             %Err = max(E1, E2);
%            [a,b,4;R,eR,abs(R*eR);N,eN,abs(N*eN);G,eG,abs(G*eG)]
         end
%     end
end 
