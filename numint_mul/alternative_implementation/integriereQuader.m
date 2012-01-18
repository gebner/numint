function[summe] = integriereQuader(f, grenzen, teilungen)
%INTEGRIEREQUADER Integriere ueber einen Hyperquader;
%
% @param[in] f          Die zu integrierende Funktion
% @param[in] grenzen    Die Grenzen des Hyperquaders:
%                       Im Feld (k, 1) ist die untere, im Feld (k,2)
%                       die obere Grenze in der k-ten Dimension
%                       gegeben;
%                       die Dimension wird an der Laenge von grenzen
%                       bestimmt.
% @param[in] teilungen  die Anzahl der bisher erfolgten Teilungen des
%                       Integrationsbereichs
% @param[out] summe     Das Integral der Funktion ueber den Quader

n = size(grenzen, 1);

summe = werteAusUndSummiere(f, grenzen, [], 6);
summe_doppelt = werteAusUndSummiere(f, grenzen, [], 8);

if log10(abs((summe_doppelt/summe)-1)) < -12 + teilungen
    summe = summe_doppelt;
    return
else
    t = mod(teilungen,n)+1;
    grenzenneu1 = grenzen;
    grenzenneu2 = grenzen;

    median = (grenzen(t,1) + grenzen(t,2))/2;
    grenzenneu1(t, 1) = median;
    grenzenneu2(t, 2) = median;
    
    summe = integriereQuader(f, grenzenneu1, teilungen+1);
    summe = summe + integriereQuader(f, grenzenneu2, teilungen+1);
end

end