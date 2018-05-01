#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
using namespace std;

const long n_max = 851;
const long infini = 1399999999;
const long iter_initialization = 1000; // Se proponen 1000 iteraciones

typedef long  type_vecteur[n_max];
typedef long type_matrice[n_max][n_max];


enum booleen {faux, vrai};


long max(long a, long b) {if (a > b) return(a); else return(b);};
double max(double a, double b) {if (a > b) return(a); else return(b);}
long min(long a, long b) {if (a < b) return(a); else return(b);}
double min(double a, double b) {if (a < b) return(a); else return(b);}
void swap(long &a, long &b) {long temp = a; a = b; b = temp;}

double temps() {return(double(clock())/double(1000));}

void a_la_ligne(ifstream & archivo){
	char poubelle[1000]; 
	archivo.getline(poubelle, sizeof(poubelle));
	}

/* Generar numeros randoms */

const long m = 2147483647; const long m2 = 2145483479; 
const long a12 = 63308; const long q12 = 33921; const long r12 = 12979; 
const long a13 = -183326; const long q13 = 11714; const long r13 = 2883; 
const long a21 = 86098; const long q21 = 24919; const long r21 = 7417; 
const long a23 = -539608; const long q23 = 3976; const long r23 = 2071;
const double invm = 4.656612873077393e-10;
long x10 = 12345, x11 = 67890, x12 = 13579, 
     x20 = 24680, x21 = 98765, x22 = 43210;

long optimo=0,menor=100000000000000;

double mon_rand()
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

long unif(long low, long high)
 {return(low + long(double(high - low + 1) *  mon_rand() ));
 }



void leer_archivo(long &n, type_matrice &a, type_matrice &b)
 {ifstream archivo;
  char nom_archivo[30];
  long i, j;

  cout << "Nombre del archivo : \n";
  cin >> nom_archivo;
 
  archivo.open(nom_archivo);
  archivo >> n; a_la_ligne(archivo);
  for (i = 1; i <= n; i = i+1) for (j = 1; j <= n; j = j+1)
    archivo >> a[i][j];
  for (i = 1; i <= n; i = i+1) for (j = 1; j <= n; j = j+1)
    archivo >> b[i][j];
  archivo.close();
 }

long calc_delta_complet2(long n, type_matrice & a, type_matrice & b,
                         type_vecteur & p, long r, long s)
 {long d;
  d = (a[r][r]-a[s][s])*(b[p[s]][p[s]]-b[p[r]][p[r]]) +
      (a[r][s]-a[s][r])*(b[p[s]][p[r]]-b[p[r]][p[s]]);
  for (long k = 1; k <= n; k = k + 1) if (k!=r && k!=s)
    d = d + (a[k][r]-a[k][s])*(b[p[k]][p[s]]-b[p[k]][p[r]]) +
            (a[r][k]-a[s][k])*(b[p[s]][p[k]]-b[p[r]][p[k]]);
  return(d);
 }

long calcule_cout(long n, type_matrice & a, type_matrice & b, type_vecteur & p)
 {long i, j;
  long c = 0;
  for (i = 1; i <= n; i = i + 1) for (j = 1; j <= n; j = j + 1)
    c = c + a[i][j] * b[p[i]][p[j]];
  return(c);
 }

void tire_solution_aleatoire(long n, type_vecteur  & p)
 {long i;
  for (i = 1; i <= n; i = i+1) p[i] = i;
  for (i = 2; i <= n; i = i+1) swap(p[i-1], p[unif(i-1, n)]);
 }

void recuit(long n, type_matrice & a, type_matrice & b,
            type_vecteur & meilleure_sol, long & meilleur_cout,
            long nb_iterations){
            
 type_vecteur p;
  long i, r, s;
  long delta;
  double cpu = temps();
  long k = n*(n-1)/2, mxfail = k, nb_fail, no_iteration;
  long dmin = infini, dmax = 0;
  double t0, tf, beta, tfound, temperature;

  for (i = 1; i <= n; i = i + 1) p[i] = meilleure_sol[i];
  long Cout = calcule_cout(n, a, b, p);
  meilleur_cout = Cout;

  for (no_iteration = 1; no_iteration <= iter_initialization;no_iteration++)
   {r = unif(1, n);
    s = unif(1, n-1);
    if (s >= r) s = s+1;

    delta = calc_delta_complet2(n,a,b,p,r,s);
    if (delta > 0)
     {dmin = min(dmin, delta); dmax = max(dmax, delta);}; 
    Cout = Cout + delta;
    swap(p[r], p[s]);
   };
  t0 = dmin + (dmax - dmin)/10.0;
  tf = dmin;
  beta = (t0 - tf)/(nb_iterations*t0*tf);

  nb_fail = 0;
  tfound = t0;
  temperature = t0;
  r = 1; s = 2;
  
  for (no_iteration = 1; no_iteration < iter_initialization - nb_iterations ; no_iteration++)
    { temperature = temperature / (1.0 + beta*temperature);

      s = s + 1;
      if (s > n)
       {r = r + 1; 
        if (r > n - 1) r = 1;
        s = r + 1;
       };

      delta = calc_delta_complet2(n,a,b,p,r,s);
      if ((delta < 0) || (mon_rand() < exp(-double(delta)/temperature)) ||
           mxfail == nb_fail)
       {Cout = Cout + delta; swap(p[r], p[s]); nb_fail = 0;}
      else nb_fail = nb_fail + 1;

      if (mxfail == nb_fail) {beta = 0; temperature = tfound;};
      if (Cout < meilleur_cout)
       {meilleur_cout = Cout;
        for (i = 1; i <= n; i = i + 1) meilleure_sol[i] = p[i];
        tfound = temperature;
        }
 
   }
   optimo= meilleur_cout;
	   if(optimo<menor){
	   	optimo= meilleur_cout;	
		   menor=optimo;
	   
	   	cout << "\nIteraciones = " << no_iteration  
	             << "  Costo = " << meilleur_cout 
	             << "  Tiempo computacional = " << temps() - cpu <<  '\n';
	      
		
	  cout << "Mejor solucion encontrada : \n";
	  for (i = 1; i <= n; i = i + 1) cout << meilleure_sol[i] << ' ';
	  cout << '\n';
	 	
		}
 }


long  n, nb_iterations, nb_res, no_res;
long Cout;
type_matrice a, b;
type_vecteur p;

main()
 {leer_archivo(n, a, b);
 
  cout << "\n Ingrese num iteraciones, num vecindarios : \n";
  cin >> nb_iterations >> nb_res;
	  for(int y=0;y<30;y++){
			  for (no_res = 1; no_res <= nb_res; no_res++)
			   {
			    tire_solution_aleatoire(n, p);
			    recuit(n,a,b,p,Cout, nb_iterations);
			   }
			   optimo=1;
			   menor=100000000000000;
			   cout<<"--------------------------------------------\n";
			   getchar();getchar();
			   
		}
}
