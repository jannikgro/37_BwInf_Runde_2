#include <stdio.h>
#include <math.h>
#define MAXDREIECKE 100
#define PI 3.1415926535897932384626
#define TOLERANZ 0.00001

#define TRUE 't'
#define FALSE 'f'
typedef char bool;

// Definition der Strukturen Punkt, Vektor und Gerade fuer die analytische Geometrie
typedef struct Punkt
{
    double x;
    double y;
    double drueber_moeglich_winkel;
    double drunter_moeglich_winkel;
    double innenwinkel;
} Punkt;

typedef struct Vektor
{
    double deltax;
    double deltay;
    double laenge;
} Vektor;

typedef struct Gerade
{
    Punkt ausgangspunkt;
    Vektor richtung;
} Gerade;

// Funktion, die den Betrag eines Vektors bildet
double vektorbetrag(Vektor v)
{
    return sqrt(v.deltax * v.deltax + v.deltay * v.deltay);
}

// Funktion, die eine vektorielle Differenz zwischen zwei Punkten bildet
Vektor vektorbestimmen(Punkt p1, Punkt p2)
{
    Vektor tmp;
    tmp.deltax = p2.x - p1.x;
    tmp.deltay = p2.y - p1.y;
    return tmp;
}

// Funktion, die eine Gerade bestimmt, auf der p1 und p2 liegen
Gerade geradebestimmen(Punkt p1, Punkt p2)
{
    Gerade gerade;
    gerade.richtung = vektorbestimmen(p1, p2);
    gerade.ausgangspunkt = p1;
    gerade.richtung.laenge = vektorbetrag(gerade.richtung);
    return gerade;
}

// modifizierte Arcustangensfunktion
double myatan(Vektor v)
{
    if (v.deltax == 0 && v.deltay == 0)
        return 10000;
    else if (v.deltax == 0)
        return (v.deltay > 0) ? 0 : PI;
    else if (v.deltay == 0)
        return (v.deltax > 0) ? PI / 2 : -PI / 2;
    else
        return atan2(v.deltax, v.deltay);
}

// Funktion, die einen Vektor umkehrt
Vektor vektorumkehren(Vektor v)
{
    v.deltax = -v.deltax;
    v.deltay = -v.deltay;
    return v;
}

// Funktion, die den Betrag einer Gleitkommazahl zurueckgibt
double betrag(double eingabe)
{
    return (eingabe > 0) ? eingabe : -eingabe;
}

// Struktur, die fuer die Beschreibung eines polygonfoermigen Hindernisses benoetigt wird
typedef struct Dreieck
{
    Punkt ecken[3];
    Gerade strecken[3];
    int modus;
} Dreieck;

int n_dreiecke = 0;            // tatsaechliche Anzahl der Dreiecke
Dreieck dreiecke[MAXDREIECKE]; // Array, in dem alle Dreiecke gespeichert werden

// Funktion, die die Informationen aus der Textdatei einliesst
void einlesen(char *dateiname)
{
    FILE *fz = fopen(dateiname, "r");
    fscanf(fz, "%d", &n_dreiecke);
    for (int i = 0; i < n_dreiecke; i++)
    {
        int test = 0;
        fscanf(fz, "\n%d", &test);
        if (test != 3)
            printf("Alle Dreiecke muessen drei Ecken haben!\n");
        for (int j = 0; j < 3; j++)
        {
            int x = 0;
            int y = 0;
            fscanf(fz, " %d %d,", &x, &y);
            dreiecke[i].ecken[j].x = x;
            dreiecke[i].ecken[j].y = y;
        }
    }
    fclose(fz);
}

// Funktion, die fuer das einzelne Dreiecke Seiten und Winkel ermittelt
void einzel_dreieck_sw(int i)
{
    for (int j = 0; j < 3; j++)
        dreiecke[i].strecken[j] = geradebestimmen(dreiecke[i].ecken[j], dreiecke[i].ecken[(j + 1) % 3]);
    for (int j = 0; j < 3; j++)
    {
        dreiecke[i].ecken[j].drueber_moeglich_winkel = myatan(dreiecke[i].strecken[j].richtung);
        dreiecke[i].ecken[j].drunter_moeglich_winkel = myatan(vektorumkehren(dreiecke[i].strecken[(j == 0) ? 2 : j - 1].richtung));
        if (dreiecke[i].ecken[j].drueber_moeglich_winkel > dreiecke[i].ecken[j].drunter_moeglich_winkel)
            dreiecke[i].ecken[j].innenwinkel = dreiecke[i].ecken[j].drueber_moeglich_winkel - dreiecke[i].ecken[j].drunter_moeglich_winkel;
        else
            dreiecke[i].ecken[j].innenwinkel = 2 * PI - betrag(dreiecke[i].ecken[j].drueber_moeglich_winkel) - betrag(dreiecke[i].ecken[j].drunter_moeglich_winkel);
    }
}

// Funktion, die ein Dreieck gegen den Uhrzeigersinn ausrichtet
void gegen_uhrzeiger(int i)
{
    int min_ecke_index = 0;
    for (int j = 1; j < 3; j++)
        if (dreiecke[i].ecken[j].y < dreiecke[i].ecken[min_ecke_index].y)
            min_ecke_index = j;
    if (dreiecke[i].ecken[min_ecke_index].drueber_moeglich_winkel < dreiecke[i].ecken[min_ecke_index].drunter_moeglich_winkel)
    {
        Dreieck ersatzdreieck;
        for (int j = 0; j < 3; j++)
            ersatzdreieck.ecken[j] = dreiecke[i].ecken[2 - j];
        dreiecke[i] = ersatzdreieck;
        einzel_dreieck_sw(i);
    }
}

// Funktion, die die Ecke mit dem kleinsten Innenwinkel ermittelt und das Dreieck so ausrichtet, dass die Ecke mit dem kleinsten Innenwinkel den Index 0 besitzt
void kleinsten_winkel_deklarieren(int i)
{
    int minindex = 0;
    for (int j = 0; j < 3; j++)
        if (dreiecke[i].ecken[j].innenwinkel < dreiecke[i].ecken[minindex].innenwinkel)
            minindex = j;
    Dreieck neuesdreieck;
    for (int j = 0; j < 3; j++)
        neuesdreieck.ecken[j] = dreiecke[i].ecken[(j + minindex) % 3];
    dreiecke[i] = neuesdreieck;
    einzel_dreieck_sw(i);
}

// Funktion, die die Dreiecke entsprechend der Vorgaben ausrichtet
void dreieckstrecken()
{
    for (int i = 0; i < n_dreiecke; i++)
    {
        einzel_dreieck_sw(i);
        gegen_uhrzeiger(i);
        kleinsten_winkel_deklarieren(i);
        dreiecke[i].modus = 0;
    }
}

short besteanordnung[MAXDREIECKE];      // In diesem Array wird die beste bisher bekannte Reihenfolge der Dreiecke gespeichert
short besteanordnung_modi[MAXDREIECKE]; // In diesem Array werden die Modi der Dreiecke der besten bisher bekannten Reihenfolge gespeichert
short laufanordnung[MAXDREIECKE];       // In diesem Array wird die aktuell ausprobierte Reihenfolge der Dreiecke gespeichert

double geringster_abstand = 10000; // Hier wird der geringste bisher erreichte Abstand gespeichert
bool used[MAXDREIECKE];            // Hier wird gespeichert, welche Dreiecke bereits fuer eine Anordnung verwendet wurden

// Hier werden die Dreiecke fuer den Beginn sinnvoll sortiert
void dreieckesortieren()
{
    for (int i = 0; i < n_dreiecke; i++)
    {
        for (int j = 0; j < n_dreiecke - 1 - i; j++)
        {
            if (dreiecke[j].ecken[0].innenwinkel > dreiecke[j + 1].ecken[0].innenwinkel)
            {
                Dreieck tmp = dreiecke[j];
                dreiecke[j] = dreiecke[j + 1];
                dreiecke[j + 1] = tmp;
            }
        }
    }
    for (int i = 0; i < n_dreiecke / 2; i += 2)
    {
        Dreieck tmp = dreiecke[i];
        dreiecke[i] = dreiecke[n_dreiecke - i - 1];
        dreiecke[n_dreiecke - i - 1] = tmp;
    }
}

// Diese Funktion positioniert ein Dreieck mit gegebenen Parametern
void positionieren(int index, double ey, double ex, double linkswinkel)
{
    double rechtswinkel = linkswinkel + dreiecke[index].ecken[0].innenwinkel;
    if (ey > TOLERANZ)
    {
        rechtswinkel = linkswinkel + PI;
        linkswinkel = rechtswinkel - dreiecke[index].ecken[0].innenwinkel;
    }
    dreiecke[index].ecken[0].x = ex;
    dreiecke[index].ecken[0].y = ey;
    dreiecke[index].ecken[1].x = ex + dreiecke[index].strecken[0].richtung.laenge * sin(rechtswinkel);
    dreiecke[index].ecken[1].y = ey + dreiecke[index].strecken[0].richtung.laenge * cos(rechtswinkel);
    dreiecke[index].ecken[2].x = ex + dreiecke[index].strecken[2].richtung.laenge * sin(linkswinkel);
    dreiecke[index].ecken[2].y = ey + dreiecke[index].strecken[2].richtung.laenge * cos(linkswinkel);
}

// Diese Funktion bestimmt den Winkel, der das vorherige Dreieck auf der rechten Seite begrenzt
double vorherrechtswinkel(int i)
{
    int preindex = laufanordnung[i - 1];
    int vorherrelevantecke = (dreiecke[preindex].ecken[0].y < TOLERANZ) ? 0 : 2;
    return myatan(vektorbestimmen(dreiecke[preindex].ecken[vorherrelevantecke], dreiecke[preindex].ecken[(vorherrelevantecke + 1) % 3]));
}

// Fuer die Anordnung mit der Spitze nach unten ermittelt diese Funktion den Winkel, unter dem das Dreieck gerade noch angeordnet werden darf (zweite Ueberschneidungspraevention)
double winkel_gegen_ueberschneidung(int j, double l, double ex)
{
    int hilfsindex = laufanordnung[j];
    int hilfsrelevant = (dreiecke[hilfsindex].ecken[0].y < TOLERANZ) ? 0 : 2;
    double qx = ex - dreiecke[hilfsindex].ecken[hilfsrelevant].x;
    if (qx == 0)
        return -PI / 2;
    double qy = -dreiecke[hilfsindex].ecken[hilfsrelevant].y;
    Vektor v = vektorbestimmen(dreiecke[hilfsindex].ecken[hilfsrelevant], dreiecke[hilfsindex].ecken[(hilfsrelevant + 1) % 3]);
    double unter_wurzel = pow(v.deltax, 2) * (-pow((qy * v.deltax - qx * v.deltay), 2) + pow(l, 2) * (pow(v.deltax, 2) + pow(v.deltay, 2)));
    double wurzel = sqrt(unter_wurzel);
    double tmpdivisor = pow(v.deltax, 2) + pow(v.deltay, 2);
    double faktor = qx + qy * v.deltax * v.deltay / tmpdivisor - qx * pow(v.deltay, 2) / tmpdivisor + wurzel / tmpdivisor;
    faktor /= v.deltax;
    if (faktor > 0 && faktor < 1)
    {
        double divident = v.deltax * (qy * v.deltax * v.deltay - qx * pow(v.deltay, 2) + wurzel);
        double divisor = qy * pow(v.deltax, 3) - v.deltay * (qx * pow(v.deltax, 2) + wurzel);
        double quotient = divident / divisor;
        double newangle = -atan(quotient);
        return newangle;
    }
    return -PI / 2;
}

// Funktion, die fuer die Anordnung mit der Spitze nach oben den Winkel ermittelt, unter dem das Dreieck gerade noch angeordnet werden darf
double konfliktwinkel(int i, int vorherrelevantecke)
{
    int index = laufanordnung[i];
    int preindex = laufanordnung[i - 1];
    double startwinkel = -PI / 2;
    for (int j = i - 2; i >= 0; i--)
    {
        for (int z = 0; z < 3; z++)
        {
            Punkt hilfspunkt = dreiecke[preindex].ecken[vorherrelevantecke];
            int hilfsindex = laufanordnung[j];
            Gerade test = geradebestimmen(hilfspunkt, dreiecke[hilfsindex].ecken[z]);
            double testwinkel = myatan(test.richtung);
            if (testwinkel > startwinkel && dreiecke[hilfsindex].ecken[z].y > TOLERANZ)
                startwinkel = myatan(vektorbestimmen(dreiecke[preindex].ecken[vorherrelevantecke], dreiecke[hilfsindex].ecken[z]));
        }
    }
    return startwinkel;
}

// Funktion, die die Summe aller verbleibenden Innenwinkel der nullten Ecke bildet
double restwinkelberechnen()
{
    double winkel = 0;
    for (int i = 0; i < n_dreiecke; i++)
        if (used[i] == FALSE)
            winkel += dreiecke[i].ecken[0].innenwinkel;
    return winkel;
}

// Funktion, die ein Dreieck mit der Spitze nach unten entsprechend der in der Dokumentation erklaerten Regeln anordnet
double spitze_unten(int i)
{
    int index = laufanordnung[i];
    dreiecke[index].modus = 3;
    double min_innenwinkel = dreiecke[index].ecken[0].innenwinkel;
    double startwinkel = -PI / 2;
    double ex = 0;
    if (i > 0)
    {
        int preindex = laufanordnung[i - 1];
        int vorherrelevantecke = (dreiecke[preindex].ecken[0].y < TOLERANZ) ? 0 : 2;
        startwinkel = vorherrechtswinkel(i);
        double restwinkel = restwinkelberechnen() + dreiecke[index].ecken[0].innenwinkel;
        if (restwinkel + startwinkel < PI / 2)
            startwinkel = (PI / 2) - restwinkel;
        ex = -dreiecke[preindex].ecken[vorherrelevantecke].y * tan(startwinkel) + dreiecke[preindex].ecken[vorherrelevantecke].x;

        // In den nachfolgenden Zeilen wird die erste Moeglichkeit der Ueberschneidung verhindert
        for (int j = i - 2; j >= 0; j--)
        {
            for (int z = 0; z < 3; z++)
            {
                Punkt hilfspunkt;
                int hilfsindex = laufanordnung[j];
                hilfspunkt.x = ex;
                hilfspunkt.y = 0;
                Gerade test = geradebestimmen(hilfspunkt, dreiecke[hilfsindex].ecken[z]);
                double testwinkel = myatan(test.richtung);
                if (test.richtung.laenge <= dreiecke[index].strecken[2].richtung.laenge && testwinkel > startwinkel && dreiecke[hilfsindex].ecken[z].y > TOLERANZ)
                {
                    startwinkel = myatan(vektorbestimmen(dreiecke[preindex].ecken[vorherrelevantecke], dreiecke[hilfsindex].ecken[z]));
                    ex = dreiecke[preindex].ecken[vorherrelevantecke].y * tan(betrag(startwinkel)) + dreiecke[preindex].ecken[vorherrelevantecke].x;
                }
            }
        }
        // Ab hier wird die zweite Moeglichkeit der Ueberschneidungen verhindert
        for (int j = i - 2; j >= 0; j--)
        {
            double newangle1 = winkel_gegen_ueberschneidung(j, dreiecke[index].strecken[2].richtung.laenge, ex);
            double newangle2 = winkel_gegen_ueberschneidung(j, dreiecke[index].strecken[0].richtung.laenge, ex) - dreiecke[index].ecken[0].innenwinkel;
            double newangle = (newangle1 > newangle2) ? newangle1 : newangle2;
            if (newangle > startwinkel)
                startwinkel = newangle;
        }
    }
    positionieren(index, 0, ex, startwinkel);
    return dreiecke[index].ecken[0].x;
}

double spitze_oben(int i)
{
    int index = laufanordnung[i];
    dreiecke[index].modus = 4;
    double min_innenwinkel = dreiecke[index].ecken[0].innenwinkel;
    int laengst_index = (dreiecke[index].strecken[0].richtung.laenge >= dreiecke[index].strecken[2].richtung.laenge) ? 0 : 2;
    int preindex = laufanordnung[i - 1];
    int vorherrelevantecke = (dreiecke[preindex].ecken[0].y < TOLERANZ) ? 0 : 2;
    double startwinkel = vorherrechtswinkel(i);
    double newangle = konfliktwinkel(i, vorherrelevantecke); // Hier wird geprueft, ob ein neuer Winkel zur Verhinderung von Ueberschneidungen gewaehlt werden muss
    if (newangle > startwinkel)
        startwinkel = newangle;
    double dy1 = dreiecke[index].strecken[0].richtung.laenge * cos(startwinkel);
    double dy2 = dreiecke[index].strecken[2].richtung.laenge * cos(startwinkel - min_innenwinkel);
    double dy = (dy1 > dy2) ? dy1 : dy2;
    double ex = dreiecke[preindex].ecken[vorherrelevantecke].x + (dy - dreiecke[preindex].ecken[vorherrelevantecke].y) * tan(startwinkel);
    positionieren(index, dy, ex, startwinkel);
    return ((dreiecke[index].ecken[1].y < TOLERANZ) ? dreiecke[index].ecken[1].x : dreiecke[index].ecken[2].x);
}

// Diese Funktion ordnet Dreiecke mit bekannten Modi an
double standardanordnung(int ziel)
{
    double strecke = 0;
    for (int i = 0; i < ziel; i++)
        used[laufanordnung[i]] = FALSE;
    for (int i = 0; i < ziel; i++)
    {
        int index = laufanordnung[i];
        used[index] = TRUE;
        int modus = dreiecke[index].modus;
        if (modus == 3)
            strecke = spitze_unten(i);
        else if (modus == 4)
            strecke = spitze_oben(i);
    }
    return strecke;
}

// Hier werden alle zum Schreiben des Ergebnisses in die Datei notwendigen Masse gespeichert
#define zoomfaktor 5
int xoffset = 300 * zoomfaktor;
int yoffset = 50 * zoomfaktor;
int hoehe = 500 * zoomfaktor;
int breite = 1250 * zoomfaktor;
int zeichenbreite = 2 * zoomfaktor;

// Diese Funktion schreibt die Umgebung und dem Weg in eine SVG-Datei
void dateischreiben(int ziel)
{
    FILE *fp = fopen("dreieckergebnis.svg", "w");
    fprintf(fp, "<svg version=\"1.1\" viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n<g transform=\"scale(1 -1)\">\n<g transform=\"translate(0 -%d)\">\n<line id=\"y\" x1=\"0\" x2=\"0\" y1=\"0\" y2=\"%d\" fill=\"none\" stroke=\"#212121\" stroke-width=\"3\"/>\n", breite, hoehe, hoehe, hoehe);
    printf("                x1|    y1;     x2|    y2;     x3|    y3;\n");
    for (int i = 0; i < ziel; i++)
    {
        int index = besteanordnung[i];
        fprintf(fp, "<polygon id=\"P%d\" points=\"", i + 1);
        for (int j = 0; j < 3; j++)
            fprintf(fp, "%d %d ", (int)round(zoomfaktor * dreiecke[index].ecken[j].x) + xoffset, (int)round(zoomfaktor * dreiecke[index].ecken[j].y) + yoffset);
        fprintf(fp, "\" fill=\"%s\" stroke=\"null\" stroke-width=\"%d\"/>\n", (i % 2 == 0) ? "ffcc99" : "#6B6B6B", zeichenbreite);
        int miny = 0;
        for (int j = 0; j < 3; j++)
            if (dreiecke[index].ecken[j].y < dreiecke[index].ecken[miny].y)
                miny = j;
        if (dreiecke[index].ecken[miny].y > TOLERANZ)
            printf("Hilfe! Das %d. Dreieck hat einen y-Abstand von %lf. Es hat den Modus %d.\n", i, dreiecke[index].ecken[miny].y, dreiecke[index].modus);
        printf("Dreieck %2d: ", i);
        for (int j = 0; j < 3; j++)
            printf("%6.1lf|%6.1lf; ", dreiecke[index].ecken[j].x, dreiecke[index].ecken[j].y);
        printf("\n");
    }
    fprintf(fp, "<line id=\"x\" x1=\"0\" x2=\"%d\" y1=\"%d\" y2=\"%d\" stroke=\"#000000\" stroke-width=\"%d\"/>\n", breite, yoffset, yoffset, zeichenbreite);
    fprintf(fp, "<line id=\"x\" x1=\"%d\" x2=\"%d\" y1=\"%d\" y2=\"%d\" stroke=\"#000000\" stroke-width=\"%d\"/>\n", xoffset, (int)(geringster_abstand * zoomfaktor) + xoffset, yoffset - 15 * zoomfaktor, yoffset - 15 * zoomfaktor, zeichenbreite);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" transform=\"scale(1 -1)\" font-family=\"sans-serif\" font-size=\"%dpx\" fill=\"black\">%.2lfm</text>\n", (int)(geringster_abstand * zoomfaktor) / 2 + xoffset, -yoffset / 4, 20 * zoomfaktor, geringster_abstand);
    fprintf(fp, "</g>\n</g>\n</svg>");
    fclose(fp);
    printf("Eine graphische Darstellung wurde in die Datei \"dreieckergebnis.svg\" geschrieben.\n");
}

// Diese Funktion rechnet die Moeglichkeiten der unterschiedlichen Anordnungen eines Pakets durch
void durchrechnen(int depth, int ausgangsdepth, int zieldepth)
{
    if (depth == zieldepth)
    {
        double abstand = standardanordnung(depth);
        geringster_abstand = abstand;
        for (int i = 0; i < n_dreiecke; i++)
        {
            besteanordnung[i] = laufanordnung[i];
            besteanordnung_modi[i] = dreiecke[laufanordnung[i]].modus;
        }
    }
    else
    {
        double restwinkel = restwinkelberechnen();
        double vorherwinkel = (depth > 0) ? vorherrechtswinkel(depth) : -PI / 2; // Dieser Winkel wurde im Abschnitt >>Loesungsidee<< als Epsilon bezeichnet
        for (int i = ausgangsdepth; i < zieldepth; i++)
        {
            if (used[i] == FALSE)
            {
                laufanordnung[depth] = i;
                used[i] = TRUE;
                double c_weg = (vorherwinkel <= 0 || restwinkel - TOLERANZ <= PI / 2 - vorherwinkel) ? spitze_unten(depth) : spitze_oben(depth);
                if (c_weg < geringster_abstand)
                    durchrechnen(depth + 1, ausgangsdepth, zieldepth);
                used[i] = FALSE;
                laufanordnung[depth] = -1;
            }
        }
    }
}

int anfangsstapel[MAXDREIECKE]; // In diesem Array werden alle Dreiecke des Anfangsstapels gespeichert

// Diese Funktion bildet den Anfangsstapel, sortiert ihn und ruft zum Schluss die Funktion >>durchrechnen<< auf
void anfangsstapeln(int cnum, double gesamtwinkel, int depth, int ausgangsdepth, int zieldepth)
{
    if (gesamtwinkel < PI / 2 && zieldepth - depth > 1)
    {
        for (int i = cnum; i < zieldepth; i++)
        {
            used[i] = TRUE;
            anfangsstapel[depth] = i;
            anfangsstapeln(i + 1, gesamtwinkel + dreiecke[i].ecken[0].innenwinkel, depth + 1, ausgangsdepth, zieldepth);
            used[i] = FALSE;
        }
    }
    else
    {
        int sortierstapel[MAXDREIECKE];
        for (int i = 0; i < depth; i++)
            sortierstapel[i] = anfangsstapel[i];
        for (int i = 0; i < depth; i++)
        {
            for (int j = 0; j < depth - i - 1; j++)
            {
                if (dreiecke[sortierstapel[j]].strecken[0].richtung.laenge < dreiecke[sortierstapel[j + 1]].strecken[0].richtung.laenge)
                {
                    int tmp = sortierstapel[j];
                    sortierstapel[j] = sortierstapel[j + 1];
                    sortierstapel[j + 1] = tmp;
                }
            }
        }
        for (int i = 0; i < depth; i++)
            laufanordnung[i] = sortierstapel[i];
        for (int i = 0; i < depth; i++)
            dreiecke[laufanordnung[i]].modus = 3;
        standardanordnung(depth);
        durchrechnen(depth, ausgangsdepth, zieldepth);
    }
}

int main(int argc, char *argv[])
{
    int dreieckepropaket = 10;
    if (argc == 1)
    {
        printf("Bitte uebergeben Sie den Dateinamen!\n");
        return 0;
    }
    if (argc == 2)
        printf("Standardmaessig wird mit %d Dreiecken pro Paket gerechnet.\n", dreieckepropaket);
    else
    {
        dreieckepropaket = atoi(argv[2]);
        if (argc == 6)
        {
            xoffset = atoi(argv[3]) * zoomfaktor;
            breite = atoi(argv[4]) * zoomfaktor;
            hoehe = atoi(argv[5]) * zoomfaktor;
        }
        else
            printf("Standardmaessig hat das Bild beim Speichern in der Datei einen xoffset von %d die Breite %d und die Hoehe %d\n", xoffset / zoomfaktor, breite / zoomfaktor, hoehe / zoomfaktor);
    }
    // Zuerst werden die zur Initialisierung notwendigen Funktionen aufgerufen
    einlesen(argv[1]);
    dreieckstrecken();
    dreieckesortieren();
    for (int i = 0; i < MAXDREIECKE; i++)
        used[i] = FALSE;
    int pakete = n_dreiecke / dreieckepropaket + ((n_dreiecke % dreieckepropaket != 0) ? 1 : 0);
    // Berechnung der einzelnen Pakete wird angeordnet
    for (int durchlauf = 0; durchlauf < pakete; durchlauf++)
    {
        geringster_abstand = 10000;
        int ausgangsdepth = durchlauf * dreieckepropaket;
        int zieldepth = ((durchlauf + 1) * dreieckepropaket > n_dreiecke) ? n_dreiecke : (durchlauf + 1) * dreieckepropaket;
        if (durchlauf == 0)
            anfangsstapeln(0, 0, 0, ausgangsdepth, zieldepth);
        else
            durchrechnen(ausgangsdepth, ausgangsdepth, zieldepth);
        for (int i = 0; i < zieldepth; i++)
        {
            laufanordnung[i] = besteanordnung[i];
            dreiecke[laufanordnung[i]].modus = besteanordnung_modi[i];
            used[i] = TRUE;
        }
        standardanordnung(zieldepth);
        printf("Jetzt wurden die Dreiecke von %d bis %d gerechnet und der bisherige Abstand ist %.2lfm.\n", ausgangsdepth, zieldepth, geringster_abstand);
    }
    dateischreiben(n_dreiecke);
    return 0;
}
