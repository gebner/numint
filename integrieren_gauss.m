function [ n,wert, fehler ] = integrieren_gauss(f,a,b, epsilon )
%INTEGRIEREN_GAUSS Summary of this function goes here
%   Detailed explanation goes here
fehler = 1;
wert_alt= 0;
n = 2;
wert= 0;
while(fehler >epsilon)
    wert = gauss_quadratur(f,n,a,b);
    fehler = abs(wert-wert_alt);
    wert_alt = wert;
    n= n+1;
    if(n>40)
        break;
    end
end

end

