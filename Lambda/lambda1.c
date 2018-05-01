/* GRASP  con Lambda-Intercambio */

#include  <stdio.h>
#include  <stdlib.h>
#include  <time.h>
#include  <sys/time.h>
#include  <sys/resource.h>
#include  <limits.h>
#include  <sys/types.h>
#include  <malloc.h>
#include  <string.h>

#define maxlista     4000
#define Max              86
#define true               1
#define false             0

#define RAND48_MULTO  (0xe66d)
#define RAND48_MULTI1   (0xdeec)
#define RAND48_MULT2    (0x0005)
#define RAND48_ADD  (0x000b)

struct listnode{
	unsigned char valor ,i , j;
};

typedef struct listnode LISTNODE;
typedef LISTNODE tlista;

typedef struct tindice{
	unsigned char i , j ;
	int valor;
};

typedef struct tindice INDICE;
typedef INDICE indices;

struct tparejas{
	char j , l;
};
typedef struct tparejas PAREJAS;
typedef PAREJAS parejas;

struct lcosto{
	unsigned char valor , i , j , k , l ;
};
typedef struct lcosto lista;
typedef lista costos;
struct tcset{
	char set [Max] , aux;
};
typedef struct tcset tset;

unsigned int seed[32] ;
long long mejor:total[32] ;
char perm_total[32][Max] ;
int iter_total[32] ;
float tproc[32] ;
long long max_inter , n , optimo;

int retval ;
unsigned char distancia[Max][Max] , flujo[Max][Max] ;  //matrices de datos

static void next (unsigned short state [] ,unsigned short mu10 ,unsigned short mul1 ,
	unsigned short mu12 ,unsigned short add)
{
	unsigned short new_state[3] ;
	unsigned long tmp ;
	
	tmp = state[0]  * mu10 + add ;
	new_state[0] = (unsigned short) (tmp & 0xffff);
	
	tmp = (tmp >> 8* sizeof (unsigned short))
		+ state[0] * mul1
		+ state[1] *mu10;
	new_state[1] = (unsigned short) ( tmp & 0xffff) ;
	
	tmp = (tmp >> 8*sizeof(unsigned short)) ;
	+ state[0] * mul2
		+ state[1] * mul1
		+ state[2] *mul0;
	new_state[2] = (unsigned short) (tmp & 0xffff) ;
	
	memcpy(state,  new_state,  3*sizeof (unsigned short)) ;
}

unsigned long mynrand48(unsigned short state[3] , unsigned short mu10 , unsigned short mul1,unsigned short mu12, unsigned short add)
{  
	
	next(state,mul0,mul1,mul2,add);
	retun(  ( (unsigned long)state [2])  *  0x8000
		+ (  ((unsigned long) state [1]) >> 1 )
		);
	
}  

unsigned long mylrand48 (unsigned short mu10 , unsigned short mul1,
	unsigned short mu12, unsigned short add, unsigned short in_s [])
{
	return (mynrand48(in_s , mu10,mu11 ,mu12 , add));
} 
void mysrand48(long seedval, unsigned short *mu10, unsigned short *mul1,
	unsigned short *mu12, unsigned short *add, unsigned short in_s [])
{
	/*Restore default multipliers and additiver. */
	*mul0 = RAND48_MULT0;
	*mul1 = RAND48_MULT1;
	*mul2 = RAND48_MULT2;
	*add = RAND48_ADD;
	/* Setup the new state */
	in_s[0] = 0x330e;
	in_s[1] = (seedval & 0xffff);
	in_s[2] = ((seedval >> 16) & 0xffff);
}
void initconta(tset *conjunto)
{
	conjunto -> aux=1;
}
/*intercambia dos numeros*/
void cambiar (char *a, char *b)
{
	int aux;
	aux = *a;
	*a = *b;
	*b = aux;
}
void printgamma(const parejas gamma[], const int total)
{
	int i;
	for (i=1;i<=n;i++)
		printf("%i,", gamma[i].1);
}
//Inicializa un conjunto
void iniset(tset *conjunto)
{
	int i;
	for (i=1; i< conjunto->aux;i++)
		conjunto->set[i]=0;
	conjunto->aux=1;
}
//inserta un elemento a un conjunto
void addset(tset *conjunto, int x)
{
	conjunto->set[conjunto->aux]=x;
	conjunto->aux++;
}
//verifica si un elemento esta en un conjunto
int inset(tset *conjunto, int x)
{
	int i;
	for (i=1;i<conjunto->aux;i++)
		if(conjunto->set[i]==x)
			return 1;
	return 0;
}
//saca un elemnto de un conjunto
void ouset(tset *conjunto, int x)
{
	int boo = false, i, j, aux[Max];
	for (i=1;i<conjunto->aux;i++)
		if (conjunto->set[i]==x)
			boo=true;
	if(boo==true)
	{for (i=1,j=1;i<conjunto->aux;i++)
		if(conjunto->set[i]!=x)
		{
			aux[j]=conjunto->set[i];
			j++;
		}
	conjunto->aux--;
	for (i=1;i<conjunto->aux;i++)
		conjunto->set[i]=aux[i];
	}
}
void inserta(tlista lista[], unsigned char matriz[] [max])
{
	int i,j,k=1;
	for (i=1;i<=n;i++)
		for (j=i+1;j<=n;j++)
		{
			lista[k].valor=matriz[i][j];
			lista[k].i=i;
			lista[k].j=j;
			k++;
		}
}
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
void ordenaij(tlista a[], int max)
{
	int j,i;
	for (j=1;j<=max;j++)
		for (i=1;i<=max;i++)
			if (i!=n)
				if(a[i].valor==a[i+1].valor)
					if(a[i].i==a[i+1].j)
					{
						if(a[i].j>a[i+1].j)
							cambiaij(&a[i],&a[i+1]);
					}
					else
						if(a[i].i>a[i+1].i)
						cambiaij(&a[i],&a[i+1]);
}
void sortldist(tlista a[],int l,int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;
	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do{
		while(a[i].valor<x)i++;
		while(x<a[j].valor)j--;
		if(i<=j)
		{
			y=a[i].valor;
			aux1=a[i].i;
			aux2=a[i].j;
			a[i].valor=a[j].valor;
			a[i].i=a[j].i;
			a[i].i=a[j].j;
			a[j].valor=y;
			a[j].i=aux1;
			a[j].j=aux2;
			i++;
			j--;
		}
	}while(i<=j);
	if(l<j)sortldist(a,l,j);
	if(i<r)sortldist(a,i,r);
}
void sortlflujo(tlista a[],int l,int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;
	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do{
		while(a[i].valor>x) i++;
		while(x>a[j].valor) j--;
		if(i<=j)
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
	if(l<j) sortlflujo(a,l,j);
	if(i<r) sortlflujo(a,i,r);
}
void sortcostos(costos a[],int l, int r)
{
	int aux1,aux2,aux3,aux4,i,j,x,y;
	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do{
		while(a[i].valor<x) i++;
		while(x<a[i].valor) i--;
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
	if(l<j) sortcostos(a,l,j);
	if(i<r) sortcostos(a,i,r);
}
/*multiplica los flujos g(int) grandes y las distancias pequenas */
void mul(costos listacostos[],const int parametro1,tlista lista1[], tlista lista2[])
{
	int i=1;
	while(i<=parametro1)
	{
		listacosto[i].valor=lista1[i].valor*lista2[i].valor;
		listacosto[i].i=lista2[i].i;
		listacosto[i].j=lista2[i].j;
		listacosto[i].k=lista2[i].i;
		listacosto[i].l=lista2[i].j;
		i++;
	}
	sortcostos(listacosto,1,parametro1);
}
void randomiza(const costos listacosto[],parejas gamma[],tset *city,tset *facility,tset *coml,tset *comj,const int parametro2,
	unsigned short mul0,unsigned short mul1, 
	unsigned short mul2,unsigned short add,unsigned short in[])
{
	int i;
	do
		{
			i=(int)mylrand48(maul0,mul1,mul2,add,in)%(parametro2)+1;
		}
	while(i<=1&&i>=parametro2+1);
	gamma[1].j=listacosto[i].i;
	gamma[1].l=listacosto[i].k;
	gamma[2].j=listacosto[i].j;
	gamma[2].l=listacosto[i].l;
	addset(facility,listacosto[i].i);
	outset(comj,listacosto[i].i);
	addset(facility,listacosto[i].j);
	outset(comj,listacosto[i].j);
	addset(city,listacosto[i].k);
	outset(coml,listacosto[i].k);
	addset(city,listacosto[i].l);
	outset(coml,listacosto[i].l);
}
/*Metodo de ordenacion quicksort*/
void sortc_gamma(indices a[],int l, int r)
{
	int aux1,aux2,i,j,x,y;
	i=l;
	j=r;
	x=a[(l+r)/2].valor;
	do
		{
			while(a[i].valor<x) i++;
			while(x<a[j].valor) j--;
			if(i<=j)
			{
				y=a[i].valor;
				aux1=a[i].j;
				aux2=a[i].i;
				a[i].valor=a[j].valor;
				a[i],j=a[j].j;
				a[i].i=a[j].i;
				a[j].valor=y;
				a[j].j=aux1;
				a[j].i=aux2;
				i++;
				j--;
			}
		}while ( i<=j);
	if(l<j) sortc_gamma(a,l,j);
	if(i<r) sortc_gamma(a,i,r);
}
/*funcion que calcula el costo de una permutacion*/
long long costo_total(const parejas g[])
{
	int i,j;
	long long suma;
	int t=0; 
	suma=0;
	for(i=1;i<=n;i++)
	{
		for(j=1;j<=n;j++)
		{
			suma=suma+(flujo[i][j])*(distancia[g[i].l][g[j].l]); 
		}
		t+=suma;
		suma=0;
	}
	return t/2;
}
/*Funcion que devuelve el numero de permutaciones que tiene el mismo costo*/
int desempate(int num, indices arr[])
{
	int i=0, j=0;
	while (num==arr[i+1].valor)
	{
		i++;
		j++;
	}
	return j;
}
/*funcion que intercambia dos registros usada por los diferentes metodos de ordenacion*/
void xchange(indices *a, indices *b)
{
	indices aux;
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
/*ordena de menor a mayor el arreglo c_gamma*/
void ordena (indices *a[maxlista])
{
	int j,i;
	for(j=1;j<=n;j++)
		for(i=1;i<=n;i++)
			if(i!=n)
				if(a[i]->valor==a[i+1]->valor)
					if(a[i]->i==a[i+1]->i)
						if(a[i]->j>a[i+1]->j)
							xchange(a[i],a[i+1]);
						else
							if(a[i]->i>a[i+1]->i)
							xchange(a[i],a[i+1]);
}
/*Construye la permutacion*/
void stage2(parejas gamma[Max],tset * city,tset *facility,tset *coml,tset *comj,
	unsigned short mul0, unsigned short mul1,unsigned short mul2,unsigned short add,unsigned short in[])
{
	int num,aux;
	long long suma;
	int i,k,j,rep;
	int m,s;
	indices c_gamma[10000];
	num=2;
	for(rep=3;rep<=n;rep++)
	{
		m=0;
		aux=0;
		for(i=1;i<=n;i++)
		{
			suma=0;
			for(k=1;k<=n;k++)
				if(inset(comj.i)&&inset(coml k))
				{
					for(j=1;j<=num;j++)
						suma=suma+((flujo[i][gamma[j].j])*distancia[k][gamma[j].l]));
						aux;
						c_gamma[aux].valor=suma;
						c_gamma[aux].i=i;
						c_gamma[aux].j=k;
						m++;
						suma=0;
				}
		}
		sortc_gamma(c_gamma,1,m);
		s=(int)mylrand48(mul0,mul1,mul2,add,in)%(desempate(c_gamma[1],valor,c_gamma))+1;
		num++;
		addset(facility,c_gamma[s].i);
		outset(comj,c_gamma[s].i);
		addset(city,c_gamma[s].j);
		outset(coml,c_gamma[s].j);
		gamma[num].j=c_gamma[s].i;
		gamma[num].l=c_gamma[s].j;
	}
} 
/*Calcula los vecinos de gamma con 2 intercambio*/
int lambda(parejas gamma[],const int costo, unsigned short mul0,unsigned short mul1,
	unsigned short mul2, unsigned short add, unsigned short in[])
{
	struct tpo{
		int costo;
		parejas gamma[Max];
	};
	typedef struct tpo TIPO;
	typedef TIPO tipoaux;
	int conta,i,j,k,l,aux,menor,error=false,total;
	parejas temp1[Max];
	tipoaux iguales[100];
	tset posibles,mov;
	aux=costo;
	conta=0;
	initset(&mov);
	initset(&posibles);
	for(i=1;i<=n;i++)
	{
		addset(&posibles,i);
	}
	do{
		j=1;
		k=j+1;
		menor=false;
		do{
			if(k==n+1)
			{
				j++; 
				k=j+1;
			}
			if ((inset(posibles.j)&&!inset(&mov,j)&&!inset(&mov,k))
				||(inset(&posibles,k)&&!inset(&mov,k)&&!inset(&mov,j)))
			{
				cambia(&gamma[j].l,&gamma[k].l);
				total=costo_total(gamma);
				if (total<aux)
				{
					conta=0;
					menor=true;
					aux=total;
					iguales[conta].costo=total;
					for(l=1;l<=n;l++)
					{
						iguales[conta].gamma[l].j=gamma[l].j;
						iguales[conta].gamma[l].j=gamma[l].l;
					}
					conta++;
				}
				else 
					if (aux==total)
					{
						iguales[conta].costo=total;
						for(l=1;l<=n;l++)
						{
							iguales[conta].gamma[l].j=gamma[l].j;
							iguales[conta].gamma[l].j=gamma[l].l;
						}
						conta++;
					}
				cambia(&gamma[k].l,&gamma[j].l);
			}
			k++;
		} while (j!=n-1);
		if (conta==0)
		{
			total=aux;
			error=true;
			break;
		}
		else
			{
				l=rand()%conta;
			}
		for(i=1;i<=n;i++)
			if (gamma[i].l!=iguales[l].gamma[i].l)
			{
				addset(&mov,i);
				outset(&posibles,i);
			}
		if (j==n-1&&k==n+1&&menor==false)
			break;
		else
			if (menor!=false)
			{
				for(i=1;i<=n;i++)
				{
					gamma[i].j=iguales[l].gamma[i].j;
					gamma[i].j=iguales[l].gamma[i].l;
				}
			}
	}while(l&&error==false);
	if (error!=true)
	{
		for (i=1;i<=n;i++)
		{
			gamma[i].j=iguales[l].gamma[i].j;
			gamma[i].j=iguales[l].gamma[i].l;
		}
		return iguales[l].costo;
	}
}

/* inicializa los conjuntos coml y comj */
void init_conjunto(tset *city, tset *facility, tset *coml, tset *comj)
{
	int i;
	initset(facility);
	initset(city);
	initset(comj);
	initset(coml);
	for (i=1;i<=n;i++)
	{
		addset(coml,i);
		addset(comj,i);
	}
}

/* Escribe los datos en un archivo de tipo de binario */
void escribe ()
{
	FILE *arch;
	int aux;
	arch=fopen("n20c.dat","wb");
	do{
		scanf("%i",&aux);
		if (aux==100)
			break;
		fwrite(&aux,sizeof(aux),1,arch);
	}while (aux!=100);
	fclose(arch);
}

/* lee de archivo binario y carga en los datos en las matrices flujo y distancia*/
void lee_archivo(const char archivo[])
{
	int i,j;
	float dato;
	FILE *arch;
		if ((arch=fopen(archivo,"r"))==NULL)
		{
			printf("error de archivo");
			exit(0);
		}
	for(i=1;i<=n;i++)
		for (j=1;j<=n;j++){
			fscanf(arch,"%d",&distancia[i][j]);
		}
	for (i=1;i<=n;i++)
		for (j=1;j<=n;j++){
			fscanf(arch,"%d",&flujo[i][j]);
		}
}

/* lee los parametros beta y alpha */ 
void lee(char str[],float *var)
{
	do
		{
			printf ("Dar el valor de %s :",str);
			scanf("%f", var);
		} while (*var<=0.0||*var>=1.0);
	printf("\n");
}

/* ordena el arreglo de menor a mayor tomando el campo 1*/
void buble(parejas gamma[])
{
	parejas aux;
	int i,j;
	for (i=1;i<=n-1;i++)
		for (j=1;j<=n-i;j++)
			if (gamma[j].j>gamma[j+1].j)
			{
				aux.j=gamma[j].j;
				aux.l=gamma[j].l;
				gamma[j].j=gamma[j+1].j;
				gamma[j].j=gamma[j+1].l;
				gamma[j+1].j=aux.j;
				gamma[j+1].j=aux.l;
			}
}

/* inicializa la premutacion poniendo 0*/
void limpiagamma(parejas gamma[])
{
	int i;
	for (i=1;i<=n;i++)
	{
		gamma[i].j=0;
		gamma[i].l=0;
	}
}

void copiap(parejas fuente[],parejas destino[])
{
	int i;
	for (i=1;i<=n;i++)
	{
		destino [i].j=fuente[i].j;
		destino [i].j=fuente[i].l;
	}
}
void grasp()
{
	static unsigned short internal_state[3] = {1, 0, 0};
	static unsigned short multiplier0 = RAND48_MULT0;
	static unsigned short multiplier1 = RAND48_MULT1;
	static unsigned short multiplier2 = RAND48_MULT2;
	static unsigned short additiver = RAND48_ADD;
	
	struct rusage ru;
	float t_inicial, t_final;
	float timep;
	time_t t;
	
	struct tcset city, facility, coml, comj, mov, posibles; //Variables para conjuntos
	parejas gamma[Max]; //Permutacion
	float beta,alpha;
	tlista lflujo[maxlista], ldist[maxlista];
	int total, res, parametro1, parametro2, band;
	char archivo[8], archivo[12] = "sal";
	costos listacosto[maxlista];
	int best, iter_local;
	unsigned long iter_grasp, iter_mejor;
	parejas mejorp[Max];
	int mejorcosto = 32000;
	int i, j, my_ide = 0;
	FILE *outfile;
	getrusage(RUSAGE_SELF, &ru);
	t_inicial = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tc_usec)/1000000);
	mysrand48(seed[my_ide], &multiplier0, &multiplier1, &multiplier2, &additiver, internal_state);
	printf("%f\n",t_inicial);
	iter_grasp = iter_mejor = 0;
	iter_local = 0;
	res = (n*(n-1))/2;
	beta = .1;
	parametro1 = (int)(parametro1*alpha)=1;
	insertar(lflujo, flujo);
	sortldist(ldist, 1, res);
	mul(listacosto, parametro1, ldist, lflujo);
	initconta(&city);
	initconta(&facility);
	initconta(&coml);
	initconta(&comj);
	
	do
		{
			limpiagamma(gamma);
			init_conjunto(&city, &facility, &coml, &comj);
			randomiza(listacosto, gamma, &city, &facility, &coml, &comj, parametro2, multiplier0, multiplier1, multiplier2, additiver, internal_state);
			stage2(gamma, &city, &facility, &coml, &comj, multiplier0, multiplier1, multiplier2, additiver, internal_state);
			buble(gamma);
			best = costo_total(gamma);
			iter_local = 0;
			do
				{
					total = lambda(gamma, best, multiplier0, multiplier1, multiplier2, additiver, internal_state);
					if(total<best)
					{	
						band=true;
						best=total;
						iter_local++;
					}
					else
						band=false;
				}while (band!=false);
			if (best<mejorcosto){
				copiap(gamma, mejorp);
				mejorcosto=best;
				iter_mejor = iter_grasp;
			}
			iter_grasp++;
		}while(((iter_grasp-iter_mejor)<max_iter)&&(optimo<mejorcosto));
	getrusage(RUSAGE_SELF, &ru);
	t_final = (float)ru.ru_utine.tv_sec + (float)(ru.ru_utine.tv_usec)/1000000;
	printf("%f\n", t_final);
	timep = (t_final - t_inicial);
	printf("%f\n", timep);
	mejor_total[my_ide] = mejorcosto;
	iter_total[my_ide] = iter_grasp;
	for(i=1; i<=n; i++){
		perm_total[my_ide][i] = mejorp[i].l;
	}
	tproc[my_ide] = timep;
	printf("posicion %d con valor %12.4f\n", my_ide, tproc[my_ide]);
	
	sprintf(archicvo, "%d", my_ide);
	strcat(archivo1, archivo);
	outfile = fopen(archivo1, "w");
	fprintf(outfile, "my_id = %7d\n", my_ide);
	for (j=1; j<=n; j++){
		fprintf(outfile, "%3d", perm_total[my_ide][j]);
	}
	fprintf (outfile, "\n");
	fprintf (outfile, "%7d", mejor_total[my_ide]);
	fprintf (outfile, "%7d", iter_total[my_ide]);
	fprintf (outfile, "%12.4f", tproc[my_ide]);
	fclose (outfile);
}


void vacia(void){
	
}


void main(argc, argv)
	int argc;
char *argv[];
{
	FILE *outfile;
	struct rusage ru;
	float t_inicial, t_final;
	float timep;
	time_t t;
	int i,j;
	
	if(argc != 7){
		printf("Forma de uso:\n");
		printf("Archivo de datos\n");
		printf("numero del problema\n");
		printf("numero maximo de iteraciones\n");
		printf("optimo\n");
		printf("numero de procesadores\n");
		printf("archivo de salida\n");
		
		exit(1213);
	}
	
	n=atoi(argv[2]);
	max_iter = atoi(argv[3]);
	optimo = atoi(argv[4]);
	Np = atoi(argv[5]);
	
	outfile = fopen(argv[6],"w");
	fprintf (outfile, "Np = %d\n", Np);
	
	for(i=0; i<Np; i++)
		seed[i]=time(&t)*(i+1);
	
	lee_archivo(argv[1]);
	
	getrusage(RUSAGE_SELF, &ru);
	t_inicial = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tv_usec)/1000000;
	
	grasp();
	
	getrusage(RUSAGE_SELF, &ru);
	t_final = (float)ru.ru_utime.tv_sec + (float)(ru.ru_utime.tv_usec)/1000000;
	
	timep = (t_final - t_final);
	
	fprintf (outfile, "Np = %d\n", Np);
	fprintf (outfile, "tiempo %12.4f\n", timep);
	
	for(i=0; i<Np; i++){
		for(j=1; j<=n; j++)
			fprintf(outfile, "%3d", perm_total[i][j]);
		fprintf (outfile, "%7d", mejor_total[i]);
		fprintf (outfile, "%7d", iter_total[i]);
		fprintf (outfile, "%12.4f", tproc[i]);
		fprintf (outfile, "\n");
	}
	fclose (outfile);
}
			






















			
			
			
			
			
