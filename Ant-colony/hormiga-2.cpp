#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <windows.h>
#include <fstream>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
const long n_max = 351;  // tamaño maximo del problema

const long infini = 1399999999;

typedef long  type_vecteur[n_max];
typedef long type_matrice[n_max][n_max];

enum booleen {faux, vrai};

void swap(long &a, long &b) {long temp = a; a = b; b = temp;}

void a_la_ligne(ifstream & archivo)
{char poubelle[1000]; archivo.getline(poubelle, sizeof(poubelle));}



/* Generar numeros randoms */
const long m = 2147483647; const long m2 = 2145483479; 
const long a12 = 63308; const long q12 = 33921; const long r12 = 12979; 
const long a13 = -183326; const long q13 = 11714; const long r13 = 2883; 
const long a21 = 86098; const long q21 = 24919; const long r21 = 7417; 
const long a23 = -539608; const long q23 = 3976; const long r23 = 2071;
const double invm = 4.656612873077393e-10;


int getMilisegundos(clock_t c)
{
    int tiempo=0;
    tiempo = (int)((c/(double)CLOCKS_PER_SEC)*1000) ;
return tiempo;

}


double rand(long & x10, long & x11, long & x12, 
            long & x20, long & x21, long & x22)
 {long h, p12, p13, p21, p23;
  h = x10/q13; p13 = -a13*(x10-h*q13)-h*r13;
  h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
  if (p13 < 0) p13 = p13 + m; if (p12 < 0) p12 = p12 + m;
  x10 = x11; x11 = x12; x12 = p12-p13; if (x12 < 0) x12 = x12 + m;
  h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
  h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
  if (p23 < 0) p23 = p23 + m2; if (p21 < 0) p21 = p21 + m2;
  x20 = x21; x21 = x22; x22 = p21-p23; if(x22 < 0) x22 = x22 + m2;
  if (x12 < x22) h = x12 - x22 + m; else h = x12 - x22;
  if (h == 0) return(1.0); else return(h*invm);
 }

long germe1 = 12345, germe2 = 67890, germe3 = 13579, 
     germe4 = 24680, germe5 = 98765, germe6 = 43210;

long unif(long low, long high)
 {return(low + long(double(high - low + 1) * 
          rand(germe1, germe2, germe3, germe4, germe5, germe6)));
 }


void leer_archivo(long & n, long & R, long & nb_iterations,  
          type_matrice & a, type_matrice & b)

 {ifstream archivo;
  char nom_archivo[30];
  long i, j;

  cout << "Nombre del archivo : \n";
  fflush(stdin);
  cin >> nom_archivo;
  archivo.open(nom_archivo);
  archivo >> n; a_la_ligne(archivo);
  cout << " num vencindarios, num iteraciones : \n";
  cin >> R >> nb_iterations; 
  if (n >= n_max) 
    {cout << "Error";

    };
  for (i = 1; i <= n; i = i+1) 
  	for (j = 1; j <= n; j = j+1)
    	archivo >> a[i][j];
  for (i = 1; i <= n; i = i+1) 
  	for (j = 1; j <= n; j = j+1)
    	archivo >> b[i][j];
  archivo.close();
 }

void imprime(long n, type_vecteur p) 
//imprime solución
 {long i; 
	 for (i = 1; i <= n; i = i + 1) 
		 cout << p[i] << ' '; cout << '\n';
 }
   
long calc_delta(long n, type_matrice & a, type_matrice & b,
                type_vecteur & p, long r, long s){
long d, k;
  d = (a[r][r]-a[s][s])*(b[p[s]][p[s]]-b[p[r]][p[r]]) +
      (a[r][s]-a[s][r])*(b[p[s]][p[r]]-b[p[r]][p[s]]);
  for (k = 1; k <= n; k = k + 1) if (k!=r && k!=s)
    d = d + (a[k][r]-a[k][s])*(b[p[k]][p[s]]-b[p[k]][p[r]]) +
            (a[r][k]-a[s][k])*(b[p[s]][p[k]]-b[p[r]][p[k]]);
  return(d);
 }

long calcule_cout(long n, type_matrice & a, type_matrice & b, type_vecteur & p)
// costo de la solución
 {long c = 0; long i, j;
  for (i = 1; i <= n; i = i + 1) for (j = 1; j <= n; j = j + 1)
    c = c + a[i][j] * b[p[i]][p[j]];
  return(c);
 }

void solution_aleatoire(long n, type_vecteur  & p)
// genera  permutación aleatoria
 {long i;
  for (i = 0; i <= n; i = i+1) p[i] = i;
  for (i = 2; i <= n; i = i+1) swap(p[i-1], p[unif(i-1, n)]);
 }

void replace(long n, type_matrice & a, type_matrice & b,
             type_vecteur & p, long & Cout)
// busqueda local
 {booleen a_tester[n_max][n_max];
  long r, s, i, ii, rr, ss, j;
  long delta;
  type_vecteur ps, pr;
  booleen ameliore = vrai;
  for (i = 1; i <= n; i = i + 1) for (j = 1; j <= n; j = j + 1)
    a_tester[i][j] = vrai;
  for (i = 1; i <= n; i = i + 1) a_tester[i][i] = faux;
  for (ii = 1; ii <= 2 && ameliore == vrai; ii = ii + 1)
   {ameliore = faux;
    solution_aleatoire(n, pr);
    for (rr = 1; rr <= n; rr = rr + 1)
     {r = pr[rr];
      solution_aleatoire(n, ps);
      for (ss = 1; ss <= n; ss = ss + 1)
       {s = ps[ss];
        if (a_tester[r][s] == vrai)
         {delta = calc_delta(n, a, b, p, r, s);
          if (delta < 0)
           {Cout = Cout + delta; swap(p[r], p[s]); ameliore = vrai;
            for (i = 1; i <= n; i = i + 1) 
              for (j = 1; j <= n; j = j + 1) a_tester[i][j] = vrai;
            for (i = 1; i <= n; i = i + 1) a_tester[i][i] = faux;
           };
          a_tester[r][s] = faux; a_tester[s][r] = faux;
         }; 
       }; 
     };
   };
 }




void initialise_trace(long n, long increment, type_matrice & trace)
// Inicialización
 {long i, j;
  for (i = 1; i <= n; i = i+1) for (j = 1; j <= n; j = j+1)
    trace[i][j] = increment;
 }

void ajourne_trace(long n, type_vecteur & p, type_vecteur & m,
                   long increment, long R, type_matrice & trace)
// Actualizar memoria
 {long i;
  for (i = 1; i <= n; i = i+1) trace[i][p[i]] = trace[i][p[i]] + increment;
  for (i = 1; i <= n; i = i+1) trace[i][m[i]] = trace[i][m[i]] + R;
 }

void tire_solution_trace(long n, type_vecteur & p, type_vecteur & m,
                         long & increment, type_matrice & trace)
// crear una solución probabilistica
 {long i, j, ii, k;
  booleen choisi[n_max];
  type_vecteur nb;
  for (i = 1; i <= n; i = i+1) nb[i] = 0;
  for (i = 1; i <= n; i = i+1) for (j = 1; j <= n; j = j+1)
    nb[i] = nb[i] + trace[i][j];
  for (i = 1; i <= n; i = i+1) choisi[i] = faux;
  long uu = unif(1, n);
  for (ii = 1; ii <= n; ii = ii+1)
   {long i = ((ii + uu) % n) + 1;
    long u = unif(1, nb[i]);
    long j = unif(1, n);
    while(choisi[j] == vrai) j = (j % n) + 1;
    long s = trace[i][j];
    while (s < u)
     {j = (j%n)+1;
      while(choisi[j] == vrai) j = (j % n) + 1;
      s = s + trace[i][j];
     };
    p[i] = j; choisi[j] = vrai;

    for (k = 1; k <= n; k = k+1) nb[k] = nb[k] - trace[k][j];
   };
  booleen identique = vrai;
  for (k = 1; k <= n && identique == vrai; k = k + 1)
    if (p[k] != m[k]) identique = faux;
  if (identique == vrai) 
   {increment = 1+increment; initialise_trace(n, increment, trace);};
 }


long  n;                                 // tamaño del problema
long Cout, meilleur_cout;                // costo de la solución actual, mejor solución
type_matrice a, b;                       // flujo y distancia matrices
type_vecteur p, mp;                      // solución actual, y mejor solución
long nb_iterations;                      // numero de iteración hormiga
type_matrice memory;                     // memoria
long increment, R;                       // parametros para el manejo de la traza
long k, no_iteration;                    // contador de iteraciones

main()
 {

 int x, contador=1;
 int t1,t2,tiempo;
 
 leer_archivo(n, R, nb_iterations, a, b);       // leer archivos y parametros
  
   
  for(x=1;x<100;x++){
	
	t1= clock();
	if (n >= n_max) return(0);

                                       // inicialización
  increment = 1;
  initialise_trace(n, increment, memory);
  meilleur_cout = infini;
                                         // iteraciones hormiga
  for (no_iteration = 1; no_iteration <= nb_iterations;
       no_iteration = no_iteration + 1){                                                  // crear nueva solución
   tire_solution_trace(n, p, mp, increment, memory);
    Cout = calcule_cout(n, a, b, p);
                                                   // busqueda local
    replace(n, a, b, p, Cout);
                                                   // checa si se mejoro
    if (Cout < meilleur_cout)
     {meilleur_cout = Cout; 
     cout<<endl;
    cout<<endl;
      cout << contador <<
        " Nuevo optimo , costo : " <<
        Cout << " Encontrado en la iteración : " << 
        no_iteration;
        contador++;
        t2=clock();
 		tiempo=getMilisegundos(t2-t1);
 		cout <<" Tiempo = "<<tiempo<< " ms"<<endl;
 		
      for (k = 1; k <= n; k = k + 1) mp[k] = p[k];
      imprime(n, p);
      increment = 1;
      initialise_trace(n, increment, memory);
     };
     
    
                                        // Actualizar
    ajourne_trace(n, p, mp, increment, R, memory);
   }
  
  getchar();getchar();
    cout<<endl<<"------------------------------------------------"<<endl; 
   }
 }

