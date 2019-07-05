#include <stdio.h>
#include <math.h>
//#include <omp.h>
#define MAXECKEN 30
#define MAXHINDERNISSE 100
#define MAXKANTE 30
#define IMPOSSIBLE -1
#define STANDARDZEIT -1000000
#define PI 3.1415926535897932384626

#define TRUE 't'
#define FALSE 'f'
typedef char bool;

// In dieser Struktur werden Zaehler und Nenner eines Bruchs fuer genaue Berechnungen gespeichert
typedef struct Bruch
{
    long long zaehler;
    long long nenner;
} Bruch;

// Euklid-Algorithmus zur Bestimmung des groessten gemeinsamen Teilers zweier Ganzzahlen
int gcd(int x, int y)
{
    return (y == 0) ? x : gcd(y, x % y);
}

// Funktion zum Kuerzen eines Bruchs
Bruch kuerzen(Bruch ausgangsbruch)
{
    if (ausgangsbruch.zaehler != 0 && ausgangsbruch.nenner != 0)
    {
        int maxkuerzen = gcd(ausgangsbruch.zaehler, ausgangsbruch.nenner);
        ausgangsbruch.zaehler /= maxkuerzen;
        ausgangsbruch.nenner /= maxkuerzen;
    }
    else
    {
        ausgangsbruch.zaehler = 0;
        ausgangsbruch.nenner = 1;
    }
    return ausgangsbruch;
}

// Funktion zur Addition zweier Brueche
Bruch addition(Bruch summand1, Bruch summand2)
{
    Bruch summe;
    summe.nenner = summand1.nenner * summand2.nenner;
    summe.zaehler = summand1.zaehler * summand2.nenner + summand2.zaehler * summand1.nenner;
    summe = kuerzen(summe);
    return summe;
}

// Multiplikation eines Bruchs mit -1
Bruch gegenbruch(Bruch ausgangsbruch)
{
    ausgangsbruch.zaehler = -ausgangsbruch.zaehler;
    return ausgangsbruch;
}

// Multiplikation zweier Brueche
Bruch multiplikation(Bruch faktor1, Bruch faktor2)
{
    Bruch produkt;
    produkt.zaehler = faktor1.zaehler * faktor2.zaehler;
    produkt.nenner = faktor1.nenner * faktor2.nenner;
    produkt = kuerzen(produkt);
    return produkt;
}

// Kehrwert eines Bruchs
Bruch umkehrung(Bruch ausgangsbruch)
{
    long long tmp = ausgangsbruch.zaehler;
    ausgangsbruch.zaehler = ausgangsbruch.nenner;
    ausgangsbruch.nenner = tmp;
    return ausgangsbruch;
}

// Konvertierung einer Ganzzahl in einen Bruch
Bruch int_to_Bruch(int ganzzahl)
{
    Bruch bruch;
    bruch.zaehler = (long long)ganzzahl;
    bruch.nenner = 1;
    return bruch;
}

// Konvertierung eines Bruchs in eine Dezimalzahl
double bruch_to_double(Bruch eingabe)
{
    return ((double)eingabe.zaehler) / ((double)eingabe.nenner);
}

// Definition der Strukturen Punkt, Vektor und Gerade fuer die analytische Geometrie
typedef struct Punkt
{
    int x;
    int y;
    bool winkelabhaengig;
    double drueber_moeglich_winkel;
    double drunter_moeglich_winkel;
} Punkt;

typedef struct Vektor
{
    int deltax;
    int deltay;
} Vektor;

typedef struct Gerade
{
    Punkt ausgangspunkt;
    Vektor richtung;
} Gerade;

// Funktion, die eine Gerade bestimmt, auf der p1 und p2 liegen
Gerade geradebestimmen(Punkt p1, Punkt p2)
{
    Gerade gerade;
    gerade.richtung.deltax = p2.x - p1.x;
    gerade.richtung.deltay = p2.y - p1.y;
    gerade.ausgangspunkt = p1;
    return gerade;
}

// Funktion, die den Betrag eines Vektors bildet
double betrag(Vektor v)
{
    return sqrt(v.deltax * v.deltax + v.deltay * v.deltay);
}

// Funktion, die prueft, ob zwei Vektoren Vielfache voneinander sind
bool vielfache(Vektor v1, Vektor v2)
{
    if (v1.deltax == 0 && v2.deltax == 0)
        return TRUE;
    else if (v1.deltay == 0 && v2.deltay == 0)
        return TRUE;
    double xverhaeltnis = ((double)(v1.deltax)) / ((double)(v2.deltax));
    double yverhaeltnis = ((double)(v1.deltay)) / ((double)(v2.deltay));
    if (xverhaeltnis == yverhaeltnis)
        return TRUE;
    else
        return FALSE;
}

double myatan(Vektor v)
{
    if (v.deltax == 0 && v.deltay == 0)
        return 10000;
    else if (v.deltax == 0)
        return (v.deltay > 0) ? 0 : PI;
    else if (v.deltay == 0)
        return (v.deltax > 0) ? PI / 2 : -PI / 2;
    else
    {
        double angle = atan2(v.deltax, v.deltay);
        if (angle > PI || angle < -PI)
            printf("WARNING:%d;%d->%lf\n", v.deltax, v.deltay, angle);
        return angle;
    }
}

unsigned long long gaussgeloest = 0;

// Funktion, die prueft, ob und wo zwei Geraden sich schneiden und damit -- wie in der Dokumentation erlaeutert -- eine direkte Verbindung unmoeglich machen
bool treffpunkt(Gerade g1, Gerade g2)
{
    gaussgeloest++;
    // Wenn die Richtungsvektoren Vielfache voneinander sind, ist die Verbindung legitim
    if (vielfache(g1.richtung, g2.richtung) == TRUE)
        return FALSE;

    Bruch var1_bruch;
    Bruch var2_bruch;
    // Abfragen zur Behandlung von speziellen Faellen, bei denen eine Vektorkomponente 0 ist
    if (g1.richtung.deltax == 0)
    {
        var2_bruch.zaehler = g1.ausgangspunkt.x - g2.ausgangspunkt.x;
        var2_bruch.nenner = g2.richtung.deltax;
        var2_bruch = kuerzen(var2_bruch);
        var1_bruch.zaehler = var2_bruch.zaehler * g2.richtung.deltay;
        var1_bruch.nenner = var2_bruch.nenner;
        var1_bruch = addition(var1_bruch, int_to_Bruch(g2.ausgangspunkt.y));
        var1_bruch = addition(var1_bruch, int_to_Bruch(-g1.ausgangspunkt.y));
        var1_bruch = multiplikation(var1_bruch, umkehrung(int_to_Bruch(g1.richtung.deltay)));
    }
    else if (g1.richtung.deltay == 0)
    {
        var2_bruch.zaehler = g1.ausgangspunkt.y - g2.ausgangspunkt.y;
        var2_bruch.nenner = g2.richtung.deltay;
        var2_bruch = kuerzen(var2_bruch);
        var1_bruch.zaehler = var2_bruch.zaehler * g2.richtung.deltax;
        var1_bruch.nenner = var2_bruch.nenner;
        var1_bruch = addition(var1_bruch, int_to_Bruch(g2.ausgangspunkt.x));
        var1_bruch = addition(var1_bruch, int_to_Bruch(-g1.ausgangspunkt.x));
        var1_bruch = multiplikation(var1_bruch, umkehrung(int_to_Bruch(g1.richtung.deltax)));
    }
    else if (g2.richtung.deltax == 0)
    {
        var1_bruch.zaehler = g2.ausgangspunkt.x - g1.ausgangspunkt.x;
        var1_bruch.nenner = g1.richtung.deltax;
        var1_bruch = kuerzen(var1_bruch);
        var2_bruch.zaehler = var1_bruch.zaehler * g1.richtung.deltay;
        var2_bruch.nenner = var1_bruch.nenner;
        var2_bruch = addition(var2_bruch, int_to_Bruch(g1.ausgangspunkt.y));
        var2_bruch = addition(var2_bruch, int_to_Bruch(-g2.ausgangspunkt.y));
        var2_bruch = multiplikation(var2_bruch, umkehrung(int_to_Bruch(g2.richtung.deltay)));
    }
    else if (g2.richtung.deltay == 0)
    {
        var1_bruch.zaehler = g2.ausgangspunkt.y - g1.ausgangspunkt.y;
        var1_bruch.nenner = g1.richtung.deltay;
        var1_bruch = kuerzen(var1_bruch);
        var2_bruch.zaehler = var1_bruch.zaehler * g1.richtung.deltax;
        var2_bruch.nenner = var1_bruch.nenner;
        var2_bruch = addition(var2_bruch, int_to_Bruch(g1.ausgangspunkt.x));
        var2_bruch = addition(var2_bruch, int_to_Bruch(-g2.ausgangspunkt.x));
        var2_bruch = multiplikation(var2_bruch, umkehrung(int_to_Bruch(g2.richtung.deltax)));
    }
    else
    {
        // Eine 3x2-Matrix wird angelegt, mit der die Loesung der Geradenschnittgleichung mithilfe eines reduzierten Gauss-Algorithmus bestimmt wird
        Bruch gaussmatrix[3][2];
        gaussmatrix[0][0] = int_to_Bruch(g1.richtung.deltax);
        gaussmatrix[0][1] = int_to_Bruch(g1.richtung.deltay);
        gaussmatrix[1][0] = int_to_Bruch(-g2.richtung.deltax);
        gaussmatrix[1][1] = int_to_Bruch(-g2.richtung.deltay);
        gaussmatrix[2][0] = int_to_Bruch(g2.ausgangspunkt.x - g1.ausgangspunkt.x);
        gaussmatrix[2][1] = int_to_Bruch(g2.ausgangspunkt.y - g1.ausgangspunkt.y);
        // Normierung der ersten Zeile darauf, dass in der ersten Zelle "1" steht
        Bruch zeilenteiler = umkehrung(gaussmatrix[0][0]);
        for (int i = 0; i < 3; i++)
            gaussmatrix[i][0] = multiplikation(gaussmatrix[i][0], zeilenteiler);
        // Elimination der ersten Variable aus der zweiten Zeile
        Bruch eliminationsfaktor = gegenbruch(gaussmatrix[0][1]);
        for (int i = 0; i < 3; i++)
            gaussmatrix[i][1] = addition(gaussmatrix[i][1], multiplikation(gaussmatrix[i][0], eliminationsfaktor));
        // Normierung der zweiten Zeile darauf, dass in der Mitte "1" steht
        zeilenteiler = umkehrung(gaussmatrix[1][1]);
        for (int i = 1; i < 3; i++)
            gaussmatrix[i][1] = multiplikation(gaussmatrix[i][1], zeilenteiler);
        // Elimination der zweiten Variable aus der ersten Zeile
        eliminationsfaktor = gegenbruch(gaussmatrix[1][0]);
        for (int i = 0; i < 3; i++)
            gaussmatrix[i][0] = addition(gaussmatrix[i][0], multiplikation(gaussmatrix[i][1], eliminationsfaktor));
        // Analyse, ob die Geraden sich so schneiden, dass ein Treffpunkt innerhalb des Bereiches vorliegt, in dem sie als Strecken tatsaechliche Bedeutung haben (Siehe Dokumentation fuer Regeln)
        var1_bruch = gaussmatrix[2][0];
        var2_bruch = gaussmatrix[2][1];
    }
    double var1 = bruch_to_double(var1_bruch);
    double var2 = bruch_to_double(var2_bruch);
    if (var1 > 0 && var1 < 1 && var2 > 0 && var2 < 1)
        return TRUE;
    else if ((var2 > 0 && var2 < 1) && (var1 == 0 || var1 == 1))
        return TRUE;
    else if (var2_bruch.zaehler == 0 && var1_bruch.zaehler < var1_bruch.nenner && var1 > 0)
        return TRUE;
    else if (var2_bruch.zaehler == var2_bruch.nenner && var1_bruch.zaehler < var1_bruch.nenner && var1 > 0)
        return TRUE;
    else
        return FALSE;
}

Vektor vektorumkehren(Vektor v)
{
    v.deltax = -v.deltax;
    v.deltay = -v.deltay;
    return v;
}

// Struktur, die fuer die Beschreibung eines polygonfoermigen Hindernisses benoetigt wird
typedef struct Hindernis
{
    int n_ecken;
    Punkt ecken[MAXECKEN];
    Gerade strecken[MAXECKEN];
} Hindernis;

// Struktur, die einen Ankunftsort auf der Strasse beschreibt (y-pos ist als Ganzzahl ein grober Naeherungswert)
typedef struct Ankunft
{
    bool possible;
    int y_pos;
    double entfernung;
} Ankunft;

// Struktur, die alle noetigen Informationen ueber einen vollstaendigen Weg vom Haus zur Strasse speichert
typedef struct Wegbeschreibung
{
    Punkt reihenfolge[MAXHINDERNISSE * MAXECKEN];
    int num_anweisungen;
    double spaeteste_zeit;
    double gesamtstrecke;
} Wegbeschreibung;

// Struktur, die die Moeglichkeiten von einem Punkt speichert oder ausschliesst
typedef struct Bitmap
{
    unsigned long long einzelspeicher[MAXHINDERNISSE * MAXECKEN / 64 + 1];
    int in_gebrauch;
} Bitmap;

typedef struct Zeitpunkt
{
    int stunden;
    int minuten;
    int sekunden;
} Zeitpunkt;

const unsigned long long A = 1; // Hilfsvariable, die zur Bitmanipulation verwendet wird

int n_hindernisse = 0;                 // tatsaechliche Anzahl der Hindernisse
Hindernis hindernisse[MAXHINDERNISSE]; // Array, in dem alle Hindernisse gespeichert werden
Punkt lisahaus;                        // Der Punkt, an dem Lisas Haus liegt

int ges_punkte;                                                            // tatsaechliche Anzahl an Ecken der Hindernisse
Punkt punkte[MAXHINDERNISSE * MAXECKEN];                                   // Array, der alle Ecken der Hindernisse speichert
double rechenmatrix[MAXHINDERNISSE * MAXECKEN][MAXHINDERNISSE * MAXECKEN]; // Speicher, in dem die Weglaengen der Wege zwischen den Ecken der Hindernisse gespeichert sind; IMPOSSIBLE bedeutet, dass der Weg nicht moeglich ist
double haus_zu_punkt[MAXHINDERNISSE * MAXECKEN];                           // Speicher, in dem die Weglaengen der Wege von Lisas Haus zu den Eckpunkten gespeichert sind; gleiche Bedeutung fuer IMPOSSIBLE
Ankunft direkt_strasse[MAXHINDERNISSE * MAXECKEN];                         // Speicher, der alle notwendigen Informationen ueber die Verbindung der Ecken zur Strasse speichert

Bitmap naechste_punkte_moeglich[MAXHINDERNISSE * MAXECKEN]; // Bitmap, die die Moeglichkeiten von einer Ecke zu anderen als "0" und "1" speichert
Bitmap von_haus_moeglich;                                   // Bitmap, die die Moeglichkeiten vom Haus zu den Ecken speichert

double min_strecke_zu_punkt[MAXHINDERNISSE * MAXECKEN]; // Array,der die kuerzeste bisher bekannte Strecke zu einem Punkt speichert

int v_lisa_km_h = 15;          // Lisas Geschwindigkeit in Kilometern pro Stunde
int v_bus_km_h = 30;           // Busgeschwindigkeit in Kilometern pro Stunde
double v_lisa;                 //Lisas Geschwindigkeit in Metern pro Sekunde
double v_bus;                  //Busgeschwindigkeit in Metern pro Sekunde
int abfahrtszeitpunkt = 27000; // Bus-Abfahrtszeitpunkt in sek (7:30 Uhr)

// Funktion, die eine leere Bitmap erstellt
Bitmap leereMap()
{
    Bitmap rueckgabe;
    rueckgabe.in_gebrauch = ges_punkte / 64 + 1;
    for (int i = 0; i < rueckgabe.in_gebrauch; i++)
        rueckgabe.einzelspeicher[i] = 0;
    return rueckgabe;
}

// Funktion, die bei einer Bitmap immer "0" und "1" tauscht
Bitmap bm_umkehren(Bitmap map)
{
    for (int i = 0; i < map.in_gebrauch; i++)
        map.einzelspeicher[i] = ~map.einzelspeicher[i];
    return map;
}

// Funktion, die den Bitoperator "&" auf zwei Bitmaps anwendet
Bitmap bit_and(Bitmap m1, Bitmap m2)
{
    for (int i = 0; i < m1.in_gebrauch; i++)
        m1.einzelspeicher[i] = m1.einzelspeicher[i] & m2.einzelspeicher[i];
    return m1;
}

// Funktion, die den Bitoperator "|" auf zwei Bitmaps anwendet
Bitmap bit_or(Bitmap m1, Bitmap m2)
{
    for (int i = 0; i < m1.in_gebrauch; i++)
        m1.einzelspeicher[i] = m1.einzelspeicher[i] | m2.einzelspeicher[i];
    return m1;
}

// Funktion, die eine Struktur vom Typ "Ankunft" initialisiert
Ankunft init_ankunft()
{
    Ankunft rueckgabe;
    rueckgabe.possible = FALSE;
    rueckgabe.y_pos = 0;
    rueckgabe.entfernung = 0;
    return rueckgabe;
}

Zeitpunkt gib_zeitpunkt(int gesamtsekunden)
{
    Zeitpunkt rueckgabe;
    rueckgabe.stunden = gesamtsekunden / 3600;
    int restsekunden = gesamtsekunden % 3600;
    rueckgabe.minuten = restsekunden / 60;
    rueckgabe.sekunden = restsekunden % 60;
    return rueckgabe;
}

// Funktion, die die Informationen aus der Textdatei einliesst
void einlesen(char *dateiname)
{
    FILE *fz = fopen(dateiname, "r");
    fscanf(fz, "%d", &n_hindernisse);
    for (int i = 0; i < n_hindernisse; i++)
    {
        hindernisse[i].n_ecken = 0;
        fscanf(fz, "\n%d", &hindernisse[i].n_ecken);
        for (int j = 0; j < hindernisse[i].n_ecken; j++)
        {
            int x = 0;
            int y = 0;
            fscanf(fz, " %d %d,", &x, &y);
            hindernisse[i].ecken[j].x = x;
            hindernisse[i].ecken[j].y = y;
        }
    }
    fscanf(fz, "\n%d %d", &lisahaus.x, &lisahaus.y);
    fclose(fz);
}

// Funktion, die fuer ein einzelnes Hindernis Strecken und Winkel bestimmt
void einzel_hindernis_sw(int i)
{
    for (int j = 0; j < hindernisse[i].n_ecken; j++)
        hindernisse[i].strecken[j] = geradebestimmen(hindernisse[i].ecken[j], hindernisse[i].ecken[(j + 1) % hindernisse[i].n_ecken]);
    for (int j = 0; j < hindernisse[i].n_ecken; j++)
    {
        hindernisse[i].ecken[j].winkelabhaengig = TRUE;
        hindernisse[i].ecken[j].drueber_moeglich_winkel = myatan(hindernisse[i].strecken[j].richtung);
        hindernisse[i].ecken[j].drunter_moeglich_winkel = myatan(vektorumkehren(hindernisse[i].strecken[(j == 0) ? (hindernisse[i].n_ecken - 1) : j - 1].richtung));
    }
}

// Funktion, die die Geradenbeschreibungen zwischen den Ecken in die Strukturen der Hindernisse eintraegt
void hindernisstrecken()
{
    //printf("\nDas sind die Hindernisse mit ihren Punkten und deren Angaben:");
    for (int i = 0; i < n_hindernisse; i++)
        einzel_hindernis_sw(i);
}

// Diese Funktion stellt sicher, dass alle Polygone im Uhrzeitersinn gedreht sind
void gegen_uhrzeiger()
{
    for (int i = 0; i < n_hindernisse; i++)
    {
        int min_ecke_index = 0;
        for (int j = 1; j < hindernisse[i].n_ecken; j++)
            if (hindernisse[i].ecken[j].y < hindernisse[i].ecken[min_ecke_index].y)
                min_ecke_index = j;
        if (hindernisse[i].ecken[min_ecke_index].drueber_moeglich_winkel < hindernisse[i].ecken[min_ecke_index].drunter_moeglich_winkel)
        {
            //printf("Hindernis %d war im Uhrzeigersinn notiert.", i + 1);
            Hindernis ersatzhindernis;
            ersatzhindernis.n_ecken = hindernisse[i].n_ecken;
            for (int j = 0; j < hindernisse[i].n_ecken; j++)
                ersatzhindernis.ecken[j] = hindernisse[i].ecken[hindernisse[i].n_ecken - j - 1];
            hindernisse[i] = ersatzhindernis;
            einzel_hindernis_sw(i);
        }
    }
}

// Funktion, die die Ecken aus den Hindernisstrukturen in den Array der Punkte eintraegt und ihre Anzahl speichert
void punkteeinfuegen()
{
    ges_punkte = 0;
    for (int i = 0; i < n_hindernisse; i++)
    {
        for (int j = 0; j < hindernisse[i].n_ecken; j++)
        {
            punkte[ges_punkte] = hindernisse[i].ecken[j];
            ges_punkte++;
        }
    }
}

// Funktion, die bestimmt, ob eine Verbindung zwischen zwei Punkten durch die Hindernisse gestoert wird
double verbindungmoeglich(Punkt p, Punkt q)
{
    Gerade weg = geradebestimmen(p, q);
    double winkel_von_startpunkt = myatan(weg.richtung);
    double winkel_von_endpunkt = myatan(vektorumkehren(weg.richtung));
    if (p.winkelabhaengig == TRUE && winkel_von_startpunkt < p.drueber_moeglich_winkel && winkel_von_startpunkt > p.drunter_moeglich_winkel && p.drunter_moeglich_winkel < p.drueber_moeglich_winkel)
        return IMPOSSIBLE;
    else if (p.winkelabhaengig == TRUE && (winkel_von_startpunkt > p.drunter_moeglich_winkel || winkel_von_startpunkt < p.drueber_moeglich_winkel) && p.drunter_moeglich_winkel > p.drueber_moeglich_winkel)
        return IMPOSSIBLE;
    else if (q.winkelabhaengig == TRUE && winkel_von_endpunkt < q.drueber_moeglich_winkel && winkel_von_endpunkt > q.drunter_moeglich_winkel && q.drunter_moeglich_winkel < q.drueber_moeglich_winkel)
        return IMPOSSIBLE;
    else if (q.winkelabhaengig == TRUE && (winkel_von_endpunkt > q.drunter_moeglich_winkel || winkel_von_endpunkt < q.drueber_moeglich_winkel) && q.drunter_moeglich_winkel > q.drueber_moeglich_winkel)
        return IMPOSSIBLE;
    else
    {
        for (int i = 0; i < n_hindernisse; i++)
            for (int j = 0; j < hindernisse[i].n_ecken; j++)
                if (treffpunkt(weg, hindernisse[i].strecken[j]) == TRUE)
                    return IMPOSSIBLE;
    }
    return betrag(weg.richtung);
}

// Funktion, die bei Uebergabe eines Punktes gemaess der 30 Grad Regel den Punkt bestimmt, bei dem Lisa bei direkter Verbindung die Strasse traefe
Punkt strasseschnitt(Punkt p)
{
    Punkt strasseerreichen;
    strasseerreichen.x = 0;
    strasseerreichen.y = p.y + (int)((double)p.x * tan(asin(v_lisa / v_bus)));
    return strasseerreichen;
}

// Funktion, die prueft, ob Lisa von ihrem Haus direkt zur Strasse gehen kann
bool lisadirekt()
{
    Punkt strasseerreichen = strasseschnitt(lisahaus);
    double haus_zu_strasse = verbindungmoeglich(lisahaus, strasseerreichen);
    if (haus_zu_strasse != IMPOSSIBLE)
    {
        double zeit = strasseerreichen.y / v_bus - haus_zu_strasse / v_lisa;
        Zeitpunkt startzeitpunkt = gib_zeitpunkt(zeit + abfahrtszeitpunkt);
        printf("Lisa kann die Strasse von ihrem Haus (%d.%d) direkt erreichen. Dazu muss sie %02d:%02d:%02d Uhr loslaufen, und die Strasse auf der y-Hoehe %d erreichen.\n", lisahaus.x, lisahaus.y, startzeitpunkt.stunden, startzeitpunkt.minuten, startzeitpunkt.sekunden, strasseerreichen.y);
        return TRUE;
    }
    else
        return FALSE;
}

// Diese Funktion fuellt saemtliche Informationen ueber die Wege zwischen den moeglichen Punkten in die dafuer vorgesehenen Arrays ein
void matrixinitialisieren()
{
    // Ab hier wird gespeichert, von welchen Punkten eine 30Grad-konforme Verbindung zur Strasse moeglich ist
    for (int i = 0; i < MAXHINDERNISSE * MAXECKEN; i++)
        direkt_strasse[i] = init_ankunft();
    for (int i = 0; i < ges_punkte; i++)
    {
        Punkt strasseerreichen = strasseschnitt(punkte[i]);
        double testentfernung = verbindungmoeglich(punkte[i], strasseerreichen);
        direkt_strasse[i].entfernung = (testentfernung != IMPOSSIBLE) ? testentfernung : sqrt(pow(strasseerreichen.y - punkte[i].y, 2) + pow(punkte[i].x, 2));
        direkt_strasse[i].y_pos = strasseerreichen.y;
        if (testentfernung != IMPOSSIBLE)
            direkt_strasse[i].possible = TRUE;
    }
    // Ab hier wird gespeichert, zu welchen Punkten Lisa von ihrem Haus direkt gehen kann; Entfernung wird gespeichert
    for (int i = 0; i < MAXHINDERNISSE * MAXECKEN; i++)
        haus_zu_punkt[i] = 0;
    for (int i = 0; i < ges_punkte; i++)
        haus_zu_punkt[i] = verbindungmoeglich(lisahaus, punkte[i]);
    // Ab hier wird gespeichert, zwischen welchen Ecken der Hindernisse eine direkte Verbindung moeglich ist; Entfernung wird gespeichert
    for (int i = 0; i < MAXHINDERNISSE * MAXECKEN; i++)
        for (int j = 0; j < MAXHINDERNISSE * MAXECKEN; j++)
            rechenmatrix[i][j] = 0;
    for (int i = 0; i < ges_punkte; i++)
    {
        for (int j = i + 1; j < ges_punkte; j++)
        {
            double entfernung = verbindungmoeglich(punkte[i], punkte[j]);
            rechenmatrix[i][j] = entfernung;
            rechenmatrix[j][i] = entfernung;
        }
    }
}

// Auf Grundlage der bereits bekannten Informationen in Rechenmatrix und Array werden hier die Bitmaps initialisiert
void bitboard_initialisieren()
{
    von_haus_moeglich = leereMap();
    for (int i = 0; i < ges_punkte; i++)
        naechste_punkte_moeglich[i] = leereMap();
    for (int i = 0; i < ges_punkte; i++)
        if (haus_zu_punkt[i] != IMPOSSIBLE)
            von_haus_moeglich.einzelspeicher[i / 64] = von_haus_moeglich.einzelspeicher[i / 64] | (A << (i % 64));
    for (int i = 0; i < ges_punkte; i++)
        for (int j = 0; j < ges_punkte; j++)
            if (rechenmatrix[i][j] != IMPOSSIBLE && i != j)
                naechste_punkte_moeglich[j].einzelspeicher[i / 64] = naechste_punkte_moeglich[j].einzelspeicher[i / 64] | (A << (i % 64));
}

// Diese Funktion druckt einen fertig ermittelten Weg zur Strasse
void printsammlung(Wegbeschreibung wb)
{
    printf("Diese %d Wege sollte Lisa von ihrem Haus {%d.%d} abgehen:\n", wb.num_anweisungen, lisahaus.x, lisahaus.y);
    for (int i = 0; i < wb.num_anweisungen; i++)
        printf("{%2d.%2d},\n", wb.reihenfolge[i].x, wb.reihenfolge[i].y);
    int gesamtzeit_bis_ankunft = wb.spaeteste_zeit + abfahrtszeitpunkt;
    Zeitpunkt startzeitpunkt = gib_zeitpunkt(gesamtzeit_bis_ankunft);
    Zeitpunkt wegdauer = gib_zeitpunkt(wb.gesamtstrecke / v_lisa);
    printf("Dann muss Lisa spaetestens um %02d:%02d:%02d Uhr loslaufen, um %.0lf Meter zurueckzulegen. Fuer den Weg benoetigt sie: %02d:%02d:%02d (Stunden:Minuten:Sekunden)\n", startzeitpunkt.stunden, startzeitpunkt.minuten, startzeitpunkt.sekunden, wb.gesamtstrecke, wegdauer.stunden, wegdauer.minuten, wegdauer.sekunden);
}

// Diese Funktion schreibt die Umgebung und dem Weg in eine SVG-Datei
void dateischreiben(Wegbeschreibung wb)
{
    FILE *fp = fopen("lisa_rennt_result.svg", "w");
    fprintf(fp, "<svg version=\"1.1\" viewBox=\"0 0 1100 750\" xmlns=\"http://www.w3.org/2000/svg\">\n<g transform=\"scale(1 -1)\">\n<g transform=\"translate(0 -750)\">\n<line id=\"y\" x1=\"0\" x2=\"0\" y1=\"0\" y2=\"750\" fill=\"none\" stroke=\"#212121\" stroke-width=\"3\"/>");
    for (int i = 0; i < n_hindernisse; i++)
    {
        fprintf(fp, "<polygon id=\"P%d\" points=\"", i + 1);
        for (int j = 0; j < hindernisse[i].n_ecken; j++)
            fprintf(fp, "%d %d ", hindernisse[i].ecken[j].x, hindernisse[i].ecken[j].y);
        fprintf(fp, "\" fill=\"#6B6B6B\" stroke=\"#212121\" stroke-width=\"2\"/>\n");
    }
    fprintf(fp, "<circle id=\"L\" cx=\"%d\" cy=\"%d\" r=\"10\" fill=\"#F42121\" stroke=\"#000080\" stroke-width=\"1\"/>\n", lisahaus.x, lisahaus.y);
    fprintf(fp, "<polyline id=\"R2\" points=\"");
    fprintf(fp, "%d %d ", lisahaus.x, lisahaus.y);
    for (int i = 0; i < wb.num_anweisungen; i++)
        fprintf(fp, "%d %d ", wb.reihenfolge[i].x, wb.reihenfolge[i].y);
    fprintf(fp, "\" fill=\"none\" stroke=\"#000080\" stroke-width=\"4\"/>\n</g>\n</g>\n</svg>");
    printf("Eine graphische Darstellung wurde in die Datei \"lisa_rennt_result.svg\" geschrieben.\n");
}

Wegbeschreibung bester_weg;             // Hier wird der beste bekannte Weg zur Strasse gespeichert
unsigned long long knotengerechnet = 0; // Hier wird die Anzahl der durchgerechneten Knoten gespeichert

// Diese Funktion uebernimmt das Ermitteln des besten Weges und speichert ihn in der oben deklarierten Struktur
void durchrechnen(Wegbeschreibung wege, Bitmap ausgeschlossene, double bisherstrecke, int index, int depth)
{
    // Hier wird die erste Optimierung implementiert: Wenn die beste theoretisch erreichbare Zeit bereits die beste bisher bekannte Zeit ueberschreitet, wird abgebrochen
    if (direkt_strasse[index].y_pos / v_bus - (bisherstrecke + direkt_strasse[index].entfernung) / v_lisa < bester_weg.spaeteste_zeit - 0.01)
        return;
    // Hier wird die dritte Optimierung implementiert: Wenn zu einem Punkt bereits eine kuerzere Strecke bekannt ist, wird abgebrochen; ansonsten wird die neue beste Zeit gespeichert
    if (min_strecke_zu_punkt[index] < bisherstrecke)
        return;
    else
        min_strecke_zu_punkt[index] = bisherstrecke;
    knotengerechnet++;                       // Die Zahl der gerechneten Knoten wird um "1" erhoeht
    wege.reihenfolge[depth] = punkte[index]; // Der aktuelle Punkt wird in die Wegbeschreibung eingetragen
    // Hier folgt die zweite Optimierung
    Bitmap auszuprobierende = bit_and(bm_umkehren(ausgeschlossene), naechste_punkte_moeglich[index]); // Hier werden die Punkte bestimmt, die in der naechsten Rechentiefe auszuprobieren sind
    ausgeschlossene = bit_or(ausgeschlossene, naechste_punkte_moeglich[index]);                       // Hier werden die Punkte bestimmt, die fuer die naechsten Rechentiefen ausgeschlossen werden
    // Wenn vom aktuellen Punkt eine direkte Verbindung zur Strasse moeglich ist, wird der neue Weg als der beste bekannte gespeichert
    if (direkt_strasse[index].possible != FALSE)
    {
        bisherstrecke += direkt_strasse[index].entfernung;
        Punkt strasseerreichen;
        strasseerreichen.x = 0;
        strasseerreichen.y = direkt_strasse[index].y_pos;
        wege.spaeteste_zeit = strasseerreichen.y / v_bus - bisherstrecke / v_lisa;
        if (wege.spaeteste_zeit > bester_weg.spaeteste_zeit)
        {
            wege.num_anweisungen = depth + 2;
            wege.reihenfolge[depth + 1] = strasseerreichen;
            wege.gesamtstrecke = bisherstrecke;
            bester_weg = wege;
        }
    }
    // Ansonsten werden durch einen rekursiven Aufruf die auszuprobierenden naechsten Punkte ausprobiert
    else
        for (int i = 0; i < ges_punkte; i++)
            if ((auszuprobierende.einzelspeicher[i / 64] & (A << (i % 64))) != 0)
                durchrechnen(wege, ausgeschlossene, bisherstrecke + rechenmatrix[index][i], i, depth + 1);
}

int main(int argc, char *argv[])
{
    v_lisa = ((double)v_lisa_km_h) / 3.6;
    v_bus = ((double)v_bus_km_h) / 3.6;
    printf("Die Geschwindigkeiten im m/s: %.3lf (Lisa), %.3lf (Bus)\n", v_lisa, v_bus);
    if (argc == 1)
    {
        printf("Bitte uebergeben Sie den Dateinamen!\n");
        return 0;
    }
    // Zuerst werden die zur Initialisierung notwendigen Funktionen aufgerufen
    einlesen(argv[1]);
    hindernisstrecken();
    gegen_uhrzeiger();
    // Ab hier wird geprueft, ob eine direkte Verbindung von Lisas Haus zur Strasse moeglich ist
    if (lisadirekt() == TRUE)
        return 0;
    // Wenn das nicht der Fall ist, werden hier die fuer die Wegsuche notwendigen Informationen gesammelt
    punkteeinfuegen();
    matrixinitialisieren();
    bitboard_initialisieren();
    // Der Speicher, der die minimale Strecke zu einem Punkt und den minimalen Gesamtweg speichert, wird so initialisiert, dass er sinnvoll verwendet werden kann
    for (int i = 0; i < ges_punkte; i++)
        min_strecke_zu_punkt[i] = -STANDARDZEIT;
    bester_weg.spaeteste_zeit = STANDARDZEIT;
    Wegbeschreibung beschreibung;
    beschreibung.spaeteste_zeit = STANDARDZEIT;
    // Von Lisas Haus werden die moeglichen Punkte ausprobiert
    for (int i = 0; i < ges_punkte; i++)
        if (haus_zu_punkt[i] != IMPOSSIBLE)
            durchrechnen(beschreibung, von_haus_moeglich, haus_zu_punkt[i], i, 0);
    // Der beste ermittelte Weg wird gedruckt
    printsammlung(bester_weg);
    dateischreiben(bester_weg);
    printf("Dafuer wurden %llu Knoten gerechnet.\n", knotengerechnet);
    printf("Dafuer wurde das Schnittgleichungssystem %llu-Male geloest.\n", gaussgeloest);
    printf("Insgesamt gibt es %d Polygomecken.\n", ges_punkte);
    return 0;
}