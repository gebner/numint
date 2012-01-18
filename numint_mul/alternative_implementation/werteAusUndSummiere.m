function[summe] = werteAusUndSummiere(f, grenzen, koordinaten, stellenzahl)
%INTEGRIEREQUADER Integriere ueber einen Hyperquader mittels Fubini
%                 und Newton-Cotes-Formeln; Reduziere dabei die
%                 Dimension des Problem rekursiv.
%
% @param[in] f           Die zu integrierende Funktion
% @param[in] grenzen     Die Integrationsgrenzen, siehe integriereQuader.m
% @param[in] koordinaten Die schon fixierten Koordinaten
% @param[in] stellenzahl Die Zahl der Stuetzstellen in einer Dimension

n = size(grenzen, 1);

a = grenzen(1,1);
b = grenzen(1,2);

koordinaten_neu = a:(b-a)/(stellenzahl-1):b;
k = length(koordinaten_neu);

w = ones(k, 1)./k;
switch k
    case 4
        w = [1/8 3/8 3/8 1/8];
    case 6
        w = [19/288 25/96 25/144 25/144 25/96 19/288];
    case 8
        w = [751/17280 3577/17280 49/640 2989/17280 ...
                2989/17280 49/640 3577/17280 751/17280];
end

w = (b-a) * w;

summe = 0;
i = 1;

if n > 1
    for x = koordinaten_neu;
        summe = summe + w(i) * werteAusUndSummiere(f, grenzen(2:end,:), ...
            [koordinaten, x], stellenzahl);
        i = i + 1;
    end
else
    for x = koordinaten_neu;
        a = f([koordinaten, x]);
        if isnan(a) || a == Inf || a == -Inf
            a = 0;
        end
        summe = summe + w(i) * a;
        i = i + 1;
    end
end

end