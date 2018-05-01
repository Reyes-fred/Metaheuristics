/* GRASP*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <task.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/prctl.h>
#include <sys/schedctl.h>
#include <string.h>

#define maxlista    3400		//maximo de elementos en listas
#define Max         83		//maximo de elementos en matriz
#define true        1
#define false       0

// constantes simbolicas para numeros seudoaleatorios

#define RAND48_MULT0		(0xe66d)
#define RAND48_MULT1		(0xdeec)
#define RAND48_MULT2		(0x0005)
#define RAND48_ADD		(0x000b)

struct listnode{
	int valor,i,j;
	};

typedef struct listnode LISTNODE;
typedef LISTNODE tlista;

typedef struct tindice{
	int i,j;
	int valor;
	};
typedef struct tindice INDICE;
typedef INDICE indices;

struct tparejas{
	   int j,l;
	   };
typedef struct tparejas PAREJAS;
typedef PAREJAS parejas;

struct lcosto{
	int valor,i,j,k,l;
	};
typedef struct lcosto lista;
typedef lista costos;

struct tcset{
	int set[Max],aux;
	};
typedef struct tcset tset;



//variables globales visibles para todos los procesadores ///////


unsigned int 		seed[32];				//semillas para cada procesador
int 			mejor_total[32];			//mejores iteraciones de cada procesador
int 			perm_total[32][Max];		//mejores permutaciones de cada procesador
int 			iter_total[32];			//mejores iteraciones para cada procesador
float 			tproc[32];			//mejores tiempos para cada procesador
int 	max_iter,						//maximo de iteraciones GRASP
	n,						//dimension de la matriz
	optimo;						//funcion objetivo
int 	Np;						//numero de procesadores
int 	retval;
int 	distancia[Max][Max],flujo[Max][Max];		//matrices de datos
int 	res,parametro1,parametro2;				//parametros para recortar las listas

/* calcula el siguiente numero seudoaleatorio*/
static void next(unsigned short state[],unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add)
{
  unsigned short new_state[3];
  unsigned long tmp;

  tmp = state[0] * mul0 + add;
  new_state[0] = (unsigned short)(tmp  & 0xffff);

  tmp = (tmp >> 8*sizeof(unsigned short))
	  + state[0] * mul1
	  + state[1] * mul0;
  new_state[1] = (unsigned short)(tmp & 0xffff);

  tmp = (tmp >> 8*sizeof(unsigned short))
	  + state[0] * mul2
	  + state[1] * mul1
	  + state[2] * mul0;
  new_state[2] = (unsigned short)(tmp & 0xffff);

  memcpy(state, new_state, 3*sizeof(unsigned short));
}

/*regresa un numero seudoaleatorio */
unsigned long mynrand48(unsigned short state[3],unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add)
{

  next(state,mul0,mul1,mul2,add);
  return( ((unsigned long)state[2]) * 0x8000
	+ ( ((unsigned long)state[1]) >> 1 )
	  );

}

/*regresa un numero seudoaleatorio */
unsigned long mylrand48(unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add,unsigned short in_s[])
{
  return(mynrand48(in_s,mul0,mul1,mul2,add));
}

/* crea una nueva secuencia de numeros seudoaleatorios*/
void mysrand48(long seedval,unsigned short *mul0,unsigned short *mul1,
				unsigned short *mul2,unsigned short *add,unsigned short in_s[])
{

  /* Restore default multipliers and additiver. */
  *mul0 = RAND48_MULT0;
  *mul1 = RAND48_MULT1;
  *mul2 = RAND48_MULT2;
  *add = RAND48_ADD;

  /* Setup the new state. */
  in_s[0] = 0x330e;
  in_s[1] = (seedval & 0xffff);
  in_s[2] = ( (seedval >> 16) & 0xffff);

}

/*intercambia dos numeros*/
void cambia(int *a,int *b)
{
	int aux;

	aux=*a;
	*a=*b;
	*b=aux;
}

/* inicializa un conjunto*/
void initset(tset *conjunto)
{
	int i;

	for (i=1;i<conjunto->aux;i++)
		conjunto->set[i]=0;
	conjunto->aux=1;
}

/*inserta un elemento a un conjunto*/
void addset(tset *conjunto,int x)
{
	conjunto->set[conjunto->aux]=x;
		conjunto->aux++;
}

/* verifica si un elemento esta en un conjunto*/
int inset(tset *conjunto,int x)
{
	int i;

	for (i=1;i<conjunto->aux;i++)
		if (conjunto->set[i]==x)
			return 1;
	return 0;
}

/* saca un elemento de un conjunto*/
void outset(tset *conjunto,int x)
{
	int boo=false,i,j,aux[Max];

	for (i=1;i<conjunto->aux;i++)
			if (conjunto->set[i]==x)
			   boo=true;
		if (boo==true)
		{
			for (i=1,j=1;i<conjunto->aux;i++)
				if (conjunto->set[i]!=x)
				{
					aux[j]=conjunto->set[i];
					j++;
				}
				conjunto->aux--;
				for (i=1;i<conjunto->aux;i++)
					conjunto->set[i]=aux[i];
		}
}

/* intercambia dos registros de lista*/
void cambiaij(tlista *a,tlista *b)
{
	  tlista aux;

	  aux.valor=a->valor;
	  aux.i=a->i;
	  aux.j=a->j;
	  a->valor=b->valor;
	  a->i=b->i;
	  a->j=b->j;
	  b->valor=aux.valor;
	  b->i=aux.i;
	  b->j=aux.j;
}

/* ordena una lista de menor a mayor tomando como referencia los indice
	i, j de la matriz */
void ordenaij(tlista a[],int max)
{
	int j,i;

	for(j=1;j<=max;j++)
		 for (i=1;i<=max;i++)
		 if (i!=max)
			if (a[i].valor==a[i+1].valor)
			   if (a[i].i==a[i+1].i)
			   {
					if (a[i].j>a[i+1].j)
						 cambiaij(&a[i],&a[i+1]);
			   }
			   else
				   if (a[i].i>a[i+1].i)
					  cambiaij(&a[i],&a[i+1]);
}

/* Ordena la lista de distancias de menor a mayor tomando como referencia
	el elemento i, j de la matriz */
void sortldist(tlista a[],int l,int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;

	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do {
		while (a[i].valor < x) i++;
		while (x < a[j].valor) j--;
		if (i<=j)
		{
			y=a[i].valor;
			aux1=a[i].i;
			aux2=a[i].j;
			a[i].valor=a[j].valor;
			a[i].i=a[j].i;
			a[i].j=a[j].j;
			a[j].valor=y;
			a[j].i=aux1;
			a[j].j=aux2;
			i++;
			j--;
		}
	}while (i<=j);
	if (l<j) sortldist(a,l,j);
	if (i<r) sortldist(a,i,r);
}

/* Ordena la lista de flujos de mayor a menor tomando como referencia
	el elemento i, j de la matriz */
void sortlflujo(tlista a[],int l,int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;

	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do {
		while (a[i].valor > x) i++;
		while (x > a[j].valor) j--;
		if (i<=j)
		{
			y=a[i].valor;
			aux1=a[i].i;
			aux2=a[i].j;
			a[i].valor=a[j].valor;
			a[i].i=a[j].i;
			a[i].j=a[j].j;
			a[j].valor=y;
			a[j].i=aux1;
			a[j].j=aux2;
			i++;
			j--;
		}
	}while (i<=j);
	if (l<j) sortlflujo(a,l,j);
	if (i<r) sortlflujo(a,i,r);
}

/* ordena la lista de costos de menor a mayor */
void sortcostos(costos a[],int l,int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;

	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do {
		while (a[i].valor < x) i++;
		while (x < a[j].valor) j--;
		if (i<=j)
		{
			y=a[i].valor;
			aux1=a[i].i;
			aux2=a[i].j;
			aux3=a[i].k;
			aux4=a[i].l;
			a[i].valor=a[j].valor;
			a[i].i=a[j].i;
			a[i].j=a[j].j;
			a[i].k=a[j].k;
			a[i].l=a[j].l;
			a[j].valor=y;
			a[j].i=aux1;
			a[j].j=aux2;
			a[j].k=aux3;
			a[j].l=aux4;
			i++;
			j--;
		}
	}while (i<=j);
	if (l<j) sortcostos(a,l,j);
	if (i<r) sortcostos(a,i,r);
}

/* calcula la lista de costos multiplicando los flujos grandes y las distancias peque¤as */
void mul_ldistlflujo(costos listacosto[],tlista ldist[],tlista lflujo[],const int parametro1)
{
	int i=1;

	while(i<=parametro1)
	{
		listacosto[i].valor=lflujo[i].valor*ldist[i].valor;
		listacosto[i].i=lflujo[i].i;
		listacosto[i].j=lflujo[i].j;
		listacosto[i].k=ldist[i].i;
		listacosto[i].l=ldist[i].j;
		i++;
	}
	sortcostos(listacosto,1,parametro1);
}

/* Toma un elemento aleatoriamente de la lista de costos
   con indices i, j, k, l y forma los primeros 2 elementos de la permutacion
   (i,k) (k,l) y los saca de sus respectivos conjuntos	*/
void randomiza(const costos listacosto[],parejas perm[],tset *city,tset *facility,tset *coml,tset *comj,const int parametro2,
				unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add,unsigned short in[])
{
	int i;

	do
	{
i=(int)mylrand48(mul0,mul1,mul2,add,in)%(parametro2)+1;	//toma un elemento //aleatoriamente
	}
	while (i<=1&&i>=parametro2+1);

	// forma los primeros dos elementos de la permutacion
	perm[1].j=listacosto[i].i;
	perm[1].l=listacosto[i].k;
	perm[2].j=listacosto[i].j;
	perm[2].l=listacosto[i].l;

	// saca los elementos asignados de los conjuntos comj y coml
	outset(comj,listacosto[i].i);
	outset(comj,listacosto[i].j);
	outset(coml,listacosto[i].k); 	outset(coml,listacosto[i].l);

	// mete los elementos asignados a los conjuntos city y facility
	addset(facility,listacosto[i].i);
	addset(facility,listacosto[i].j);
	addset(city,listacosto[i].k);
	addset(city,listacosto[i].l);
}

/*funcion que calcula el costo de una permutacion*/
int costo_total(const parejas g[])
{
	int i,j,suma;
	int t=0;
	suma=0;
	for (i=1;i<=n;i++)
	{
		for (j=1;j<=n;j++)
		{
			suma=suma+(flujo[i][j])*(distancia[g[i].l][g[j].l]);
		}
		t+=suma;
		suma=0;
	}
	return t / 2;
}

/* Funcion que devuelve el numero de elementos que tienen el mismo costo*/
int desempate(int num, indices arr[])
{
	int i=0,j=0;
	while (num==arr[i+1].valor)
	{
		i++;
		j++;
	}
	return j;
}

/* Construye la permutacion*/
void stage2(parejas perm[Max], tset *city,tset *facility,tset *coml,tset *comj,
			unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add,unsigned short in[])

{
	int num,aux,suma;
	int i,k,j,rep;
	int m,s;
	indices gamma[250000];
	num=2;
	for (rep=3;rep<=n;rep++)
	{
		m=0;
		aux=0;
		for (i=1;i<=n;i++)
		{
			suma=0;
			for (k=1; k<=n;k++)
				if (inset(comj,i)&&inset(coml,k))
				{
					for (j=1;j<=num;j++)
						suma=suma+((flujo[i][perm[j].j])*(distancia[k][perm[j].l]));
					aux++;
					gamma[aux].valor=suma;
					gamma[aux].i=i;
					gamma[aux].j=k;
					m++;
					suma=0;
				}
		}

		//ordena el arreglo gamma
		sortldist(gamma,1,m);

		//toma un elemento aleatoriamente
		s=(int)mylrand48(mul0,mul1,mul2,add,in)%(desempate(gamma[1].valor,gamma))+1;
		num++;

		//asigna el elemento a la permutacion
		perm[num].j=gamma[s].i;
		perm[num].l=gamma[s].j;

		//inserta el elemento asignado a los conjuntos facility y city
		addset(facility,gamma[s].i);
		addset(city,gamma[s].j);

		//saca el elemento asignado de los conjuntos coml y comj
		outset(coml,gamma[s].j);
		outset(comj,gamma[s].i);
	}
}

/* calcula los vecinos de la permutacion y sus costos y
devuelve la mejor permutacion y el mejor costo encontrado.
Si algun vecino tiene un costo igual al mejor encontrado, se decide el
movimiento de cambio con un volado*/
int two_exchange(parejas perm[], const int costo,
				 unsigned short mul0,unsigned short mul1,
				unsigned short mul2,unsigned short add,unsigned short in[])
{
 int moneda;
 int i,j;
 int j0,i0;
 int total,aux=costo,cambio=false;
 for (i=1;i<n;i++)
 {
	for (j=i+1;j<=n;j++)
	{
		cambia(&perm[i].l,&perm[j].l);
		total=costo_total(perm);
		moneda=1;
		if (total==optimo)
		{
			aux=total;
			return aux;
		}
		if (total==aux)

		   //volado
		   moneda=(int)mylrand48(mul0,mul1,mul2,add,in) % 2;
		if (total<=aux&&moneda==1)
		{
		   j0=j;
		   i0=i;
		   aux=total;
		   cambio=true;
		}
		cambia(&perm[j].l,&perm[i].l);
	}
 }
 if (cambio==true)
 {
	cambia(&perm[i0].l,&perm[j0].l);
	return aux;
 }
 else
	 {
		return aux;
	 }
}

/*inicializa los conjuntos coml, comj, city, facility */
void init_conjunto(tset *city,tset *facility,tset *coml,tset *comj)
{
	int i;
	initset(facility);
	initset(city);
	initset(comj);
	initset(coml);
	//inicializa los conjuntos coml y comj de 1..n
	for (i=1;i<=n;i++)
	{
		addset(coml,i);
		addset(comj,i);
	}
}

/* lee de archivo binario y carga en los datos en las matrices flujo y distancia*/
void carga_archivo(const char archivo[])
{
	int dato,i,j;
	FILE *arch;
	if ((arch=fopen(archivo,"rb"))==NULL)
	{
		printf("Error de archivo");
		exit(0);
	}
	i=j=1;
	while ((i<=n)&&(feof(arch)==0))
	{
		if (j==n+1)
		{
			i++;
			j=i;
		}
		if (i==j)
		{
			distancia[i][j]=0;
			j++;
		}
		else
		{
			fread(&dato,sizeof(dato),1,arch);
			distancia[i][j]=dato;
			if (distancia[i][j]==distancia[j][i])
				distancia[i][j]=distancia[j][i];
			distancia[j][i]=distancia[i][j];
			j++;
		}
	}
	i=j=1;
	while ((i<=n)&&(feof(arch)==0))
	{
		if (j==n+1)
		{
			i++;
			j=i;
		}
		if (i==j)
		{
			flujo[i][j]=0;
			j++;
		}
		else
		{
			fread(&dato,sizeof(dato),1,arch);
			flujo[i][j]=dato;
			if (flujo[i][j]==flujo[j][i])
				flujo[i][j]=flujo[j][i];
			flujo[j][i]=flujo[i][j];
			j++;
		}
	}
}

/* ordena la permutacion */
void sort_perm(parejas perm[])
{
	parejas aux;
	int i,j;
	for (i=1;i<=n-1;i++)
		for (j=1;j<=n-i;j++)
			if (perm[j].j>perm[j+1].j)
			{
				aux.j=perm[j].j;
				aux.l=perm[j].l;
				perm[j].j=perm[j+1].j;
				perm[j].l=perm[j+1].l;
				perm[j+1].j=aux.j;
				perm[j+1].l=aux.l;
			}
}

/* copia una permutacion de fuente a destino*/
void copiap(parejas fuente[],parejas destino[])
{ 
	int i;
	for (i=1;i<=n;i++)
	{
		destino[i].j=fuente[i].j;
				destino[i].l=fuente[i].l;
	}
}

/* carga las listas de distancias y flujos de la matriz de datos */
void carga_listas(tlista ldist[],tlista lflujo[])
{
   int i,j,k=1;
   for (i=1;i<=n;i++)
	   for(j=i+1;j<=n;j++)
	   {
		  lflujo[k].valor=flujo[i][j];
		  lflujo[k].i=i;
		  lflujo[k].j=j;
		  ldist[k].valor=distancia[i][j];
		  ldist[k].i=i;
		  ldist[k].j=j;
		  k++;
	   }
}

void grasp()
{
	//variables para funciones aleatorias
	static unsigned short internal_state[3] = {1, 0, 0};
	static unsigned short multiplier0 = RAND48_MULT0;
	static unsigned short multiplier1 = RAND48_MULT1;
	static unsigned short multiplier2 = RAND48_MULT2;
	static unsigned short additiver = RAND48_ADD;

	struct rusage ru;
	float t_inicial,t_final,timep;
	time_t t;
	struct tcset city,facility,coml,comj,mov,posibles;		//variables para conjuntos
	parejas perm[Max];  				//permutacion
	char archivo[8],archivo1[12]="sal";
	unsigned long  iter_grasp,iter_mejor,			//contadores de iteraciones
			mejorcosto = 320000;		// mejor costo
	parejas mejorp[Max];				//mejor permutacion
	int i,j,
		best,					//mejor costo de permutacion
		total,					//costo de permutacion
		band,
		my_ide = m_get_myid();			//identificador de procesador
	tlista lflujo[maxlista],ldist[maxlista];			//listas de flujos y distancias
	costos listacosto[maxlista];				//lista de costos
	FILE *outfile;					//archivo de salida de procesador

	//tomar el tiempo inicial de ejecucion correspondiente al procesador
	getrusage(RUSAGE_SELF,&ru);
	t_inicial = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tv_usec)/1000000;

	//crear una secuencia de numeros seudoaleatorios con la semilla correspondiente al procesador
	mysrand48(seed[my_ide],&multiplier0,&multiplier1,&multiplier2,&additiver,internal_state);

	//cargar la lista de distancias y flujos de la matriz de datos
	carga_listas(ldist,lflujo);

	//ordenar los flujos de mayor a menor
	sortlflujo(lflujo,1,res);
	//ordenar los flujos con respecto a los indices
	ordenaij(lflujo,res);
	//ordenar las distancias de menor a mayor
	sortldist(ldist,1,res);
	//ordenar los flujos con respecto a los indices
	ordenaij(ldist,res);
	//calcular la lista de costos
	mul_ldistlflujo(listacosto,ldist,lflujo,parametro1);
	//inicializar los contadores de iteraciones
	iter_grasp=iter_mejor=0;

	//inicializar los conjuntos
	initset(&city);
	initset(&facility); 	initset(&coml);
	initset(&comj);

	do
	{
		//inicializar con los valores iniciales a los conjuntos
		init_conjunto(&city,&facility,&coml,&comj);

		//tomar un elemento de la lista de costos aleatoriamente y crear
		// los primeros dos elementos de la permutacion
		randomiza(listacosto,perm,&city,&facility,&coml,&comj,parametro2,multiplier0,multiplier1,multiplier2,
			additiver,internal_state);
		//construir el resto de la permutacion
		stage2(perm,&city,&facility,&coml,&comj,multiplier0,multiplier1,multiplier2,additiver,internal_state);

		//ordenar la permutacion
		sort_perm(perm);

		//calcular el costo de la permutacion inicial
		best=costo_total(perm);
		do
		{
			//calcular los vecinos y encontrar el menor de ellos
			total=two_exchange(perm, best,multiplier0,multiplier1,multiplier2,additiver,internal_state);
			if (total<best)
			{
				band=true;
				best=total;
			}
			else
				band=false;
		}while (band!=false);//repetir hasta no encontrar un vecino menor
		//¨el costo del vecino es menor que el mejor encontrado
		if (best<mejorcosto)
		{
			//guardar la encontrada permutacion con costo menor
			copiap(perm,mejorp);
			//actualizar mejor costo encontrado
			mejorcosto=best;
			//actualizar mejor iteracion
			iter_mejor = iter_grasp;
		}
		iter_grasp++;
	//repetir hasta no encontrar un costo menor en (iter_grasp - iter_mejor)
	//iteraciones > max_iter o hasta encontrar el optimo
	}while(((iter_grasp-itermejor)<max_iter)&&(optimo<mejorcosto));
	//tomar el tiempo final de ejecucion correspondiente al procesador
	getrusage(RUSAGE_SELF,&ru);
	t_final = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tv_usec)/1000000;
	//calcular el tiempo total de ejecucion
	timep = (t_final- t_inicial);
	//guardar  mejor costo correspondiente al numero procesador
	mejor_total[my_ide] = mejorcosto;
	//guardar mejor tiempo correspondiente al numero procesador
	tproc[my_ide]=timep;
	//guardar numero de iteraciones correspondiente al numero procesador
	iter_total[my_ide] = iter_grasp;
	//guardar la mejor permutacion correspondiente al numero procesador
	for (i=1;i<=n;i++)
	{
		perm_total[my_ide][i] =  mejorp[i].l;
	}

	//imprime cada uno de los datos al archivo de salida de procesador
	sprintf(archivo,"%d",my_ide);
	strcat(archivo1,archivo);
	outfile = fopen (archivo1,"w");
	fprintf (outfile,"my_id = %7d\n",my_ide);
	for (j=1;j<=n;j++)
	{
		fprintf (outfile,"%3d",perm_total[my_ide][j]);
	}
	fprintf (outfile,"\n");
	fprintf (outfile,"%7d",mejor_total[my_ide]);
	fprintf (outfile,"%7d",iter_total[my_ide]);
	fprintf (outfile,"%12.8f",tproc[my_ide]);
	fclose (outfile);
}

//subrutina vacia para crear hilos
void vacia(void)
{
	   /*vacia routine to create threads */
}

void main(argc, argv)
int argc;
char *argv[];
{
	FILE *outfile;					//archivo de salida general
	struct rusage ru;
	float t_inicial,t_final,timep;
	time_t t;
	int i,j;

	if (argc != 7)
	{
		exit (1213);
	}

	carga_archivo(argv[1]);			//nombre del archivo de datos
	n = atoi(argv[2]);				//dimension de la matriz	
	max_iter = atoi(argv[3]);			//maximo de iteraciones
	optimo = atoi(argv[4]);			//funcion objetivo
	Np = atoi(argv[5]);				//numero de procesadores
	outfile = fopen (argv[6],"w");		//archivo de salida general


	m_set_procs(Np);				//fijar el numero de proceasores a utilizar
	fprintf (outfile,"Np = %d\n",Np);

	// crear semillas de numeros seudoaleatorios para cada uno de los procesadores
	for (i=0;i<Np;i++)
		seed[i] = time(&t)*(i+1);

	//calcular parametros para cortar las listas de distancias y flujos
	res=(n*(n-1))/2;
	parametro1=(int)(res*.1)+1;
	parametro2=(int)(parametro1*.5)+1;
	//tomar el tiempo inicial de ejecucion
	getrusage(RUSAGE_SELF,&ru);
	t_inicial = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tv_usec)/1000000;

	if (Np > 1)
	{
		m_fork(vacia); 		//crear un proceso inutil en cada procesador
		m_park_procs();		//suspende procesos hijos

		retval = schedctl (SCHEDMODE,SGS_GANG);
		if (retval == -1 )
		{
			perror("schedctl");
			exit(100);
		}

		m_rele_procs();		//reiniciar procesos
	}

	m_fork(grasp);		//ejecutar grasp en cada procesador
	m_kill_procs();		//matar los procesos creados

	//tomar el tiempo final de ejecucion
	getrusage(RUSAGE_SELF,&ru);
	t_final = (float)ru.ru_utime.tv_sec +(float)(ru.ru_utime.tv_usec)/1000000;

	//calcular el tiempo total de ejecucion
	timep = (t_final- t_inicial);

	//Guardar los datos en el archivo de salida general
	fprintf (outfile,"Np = %d\n",Np);
	fprintf (outfile,"tiempo %12.8f\n",timep);
	for (i=0;i<Np;i++)
	{
		for (j=1;j<=n;j++)
			fprintf (outfile,"%6d",perm_total[i][j]);
		fprintf (outfile,"%7d",mejor_total[i]);
		fprintf (outfile,"%7d",iter_total[i]);
		fprintf (outfile,"\n%12.8f",tproc[i]);
		fprintf (outfile,"\n");
	}
	fclose (outfile);
}

