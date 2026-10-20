# Wärmeleitung auf einem Stab

`index.html` direkt im Browser öffnen. Das Applet funktioniert offline ohne Build. Alternativ: `npm start`.

## Ablauf

1. Randbedingungen wählen: feste Endtemperatur, isolierte Enden oder Ringstab.
2. Anfangsfunktion eingeben oder ein passendes Beispiel wählen. Nicht erfüllte Randbedingungen blockieren den nächsten Schritt.
3. Mit dem Schieberegler 1–10 Fourierfrequenzen wählen. Gezeigt werden die Anfangsfunktion und die **Partialsumme**, keine einzelnen Moden. Direkt unter dem Regler erscheint die aktuelle Partialsumme als LaTeX-Formel.
4. Unter dem Graphen zeigt ein zylindrischer Stab die Anfangstemperatur: blau ist kalt, rot ist heiss.
5. **Zur Zeitanimation** wechseln. Beide Kurven entwickeln sich mit der Wärmegleichung: Referenzlösung und Partialsumme mit der gewählten Fourierordnung. Der Stab folgt ausschliesslich der Referenzlösung.

`alpha = 1` und `L = pi` sind fest. Alle Probleme werden auf `[-L,L]` gelöst; die Achse ist mit `−L`, `−L/2`, `0`, `L/2`, `L` beschriftet. Die Eingabe darf `L` enthalten, z. B. `1-(x/L)^2`.

Die Skala wird einmal pro Anfangsfunktion bestimmt und bleibt bei Änderungen der Fourierordnung und der Zeit konstant. Sie umfasst die Anfangsfunktion, alle auswählbaren Partialsummen und die Referenzreihe bei Zeit null. Die ursprüngliche Anfangsfunktion wird bei Zeit null direkt ausgewertet; für positive Zeiten nutzt die Referenz 200 Fourierfrequenzen. Diese ist eine numerische Referenz, keine garantierte exakte Lösung für beliebige Eingaben. Die Referenz hängt nicht vom Fourierregler ab.

Eine Änderung der Randbedingungen oder Anfangsfunktion setzt die folgenden Schritte zurück. Beim Verschieben des Fourierreglers wird eine laufende Animation pausiert; die aktuelle Zeit bleibt erhalten. Bei Neumann und periodischen Randbedingungen kommt der konstante Mittelwert zur Summe der gewählten Frequenzen hinzu. Beim Ring enthält eine Frequenzmode einen Sinus- und einen Kosinusanteil.

Die Animation dauert 30 Sekunden. Der Zeitbereich endet bei 1,5 Abklingzeiten einer anhand der nichtkonstanten spektralen Energie bestimmten Frequenz. Reine hohe Frequenzen bekommen dadurch einen entsprechend kürzeren Zeitbereich. Die Zuordnung `t = t_max * s²` für Schieberposition bzw. Animationsfortschritt `s` gibt den frühen Zeiten mehr Raum. Das Zeitfeld zeigt weiterhin die physikalische Zeit; die PDE bleibt unverändert.

Die Farbskala des Stabs wird einmal aus der Anfangsfunktion und Temperatur null bestimmt. Sie wird während der Animation und beim Ändern der Fourierordnung nicht angepasst. Bei Zeit null werden die Eingabewerte direkt verwendet, danach dieselben Referenzwerte wie für den Graphen. Zurücksetzen stellt die Anfangsfarben wieder her. Der periodische Stab wird abgerollt dargestellt.

Die Formelausgabe rundet Koeffizienten auf zwei Dezimalstellen; auf null gerundete Terme werden ausgeblendet. Berechnung, Graphen und Temperaturstab verwenden weiterhin die ungerundeten Werte. `series.js` erkennt zusätzlich endliche lineare Kombinationen der passenden Eigenfunktionen strukturell. Auch erkannte exakte Koeffizienten werden für die Anzeige auf zwei Dezimalstellen gerundet. Ein Gleichheitszeichen wird nur verwendet, wenn alle enthaltenen Frequenzen ausgewählt sind und die Rundung keinen Koeffizienten verändert; sonst steht ein Näherungszeichen. Ein kleiner numerischer Rest oder übereinstimmende Stichproben allein begründen keine exakte Gleichheit; nicht erkannte Identitäten bleiben konservativ als Näherung markiert.

## Mathematik und Prüfung

Grundlage sind die analytischen Fourierlösungen aus **Analysis_III_BAUG_V2.pdf**, Kapitel 3.6–3.7, S. 96–98. Für Dirichlet und Neumann ist das Skriptintervall auf `[-L,L]` verschoben und auf die Länge `2L` angepasst:

- Dirichlet: `sin(n*pi*(x+L)/(2*L))` mit Zeitfaktor `exp(-(n*pi/(2*L))^2*t)`.
- Neumann: `cos(n*pi*(x+L)/(2*L))` mit demselben Zeitfaktor und konstantem Mittelwert.
- Periodisch: `sin(n*pi*x/L)`, `cos(n*pi*x/L)` mit Zeitfaktor `exp(-(n*pi/L)^2*t)` und konstantem Mittelwert.

Die Koeffizienten werden mit zusammengesetzter Gauss-Legendre-Quadratur und Verfeinerungsvergleich berechnet. Kein numerisches Zeitschrittverfahren. Der Rechenkern bietet mit `centered: true` das symmetrische Intervall; ohne diese Option bleiben die ursprünglichen Intervalle erhalten.

Die Randprüfung verwendet Funktionswerte und einseitige numerische Ableitungen mit relativer Toleranz. Sie ist keine symbolische Garantie für beliebige Eingaben. Der Parser akzeptiert z. B. `sin(pi*x/L)`, `1-(x/L)^2` oder `x < 0 ? 1 : 0` und lässt keine Zuweisungen oder Objektzugriffe zu.

`npm test` prüft Fourierkoeffizienten, Randverträglichkeit, Mittelwerterhaltung, Zeitentwicklung, die PDE, die angepassten Eigenfunktionen auf `[-L,L]` und Eingabevalidierung. Lokale Bibliotheken in `vendor/`: KaTeX 0.16.11 und mathjs 14.8.1 einschließlich Lizenzen.
