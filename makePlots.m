
loglog(1./y,y,'black')
hold on
loglog(1./y,E,'blue')
xlabel 'Fehlerschranke (1/epsilon)'
loglog(1./y,F,'green')
legend('Fehlerschranke','Fehlersch�tzer','Tats�chlicher Fehler')
hold off