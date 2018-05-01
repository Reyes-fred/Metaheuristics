/* Grasp con desempate en

		stage2
		two exchange
				*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
//#include <sys/resource.h>
//#include <task.h>
#include <limits.h>
#include <sys/types.h>
//#include <sys/prctl.h>
//#include <sys/schedctl.h>
#include <malloc.h>
#include <string.h>


#define maxlista    4000
#define Max         86
#define true        1
#define false       0

#define RAND48_MULT0 (0xe66d)
#define RAND48_MULT1 (0xdeec)
#define RAND48_MULT2 (0x0005)
#define RAND48_ADD (0x000b)

   struct listnode{
      unsigned char valor,i,j;
   };

   typedef struct listnode LISTNODE;
   typedef LISTNODE tlista;

   typedef struct tindice{
      unsigned char i,j;
      int valor;
   };

   typedef struct tindice INDICE;
   typedef INDICE indices;



   struct tparejas{
      char j,l;
   };
   typedef struct tparejas PAREJAS;
   typedef PAREJAS parejas;

   struct lcosto{
      unsigned char valor,i,j,k,l;
   };
   typedef struct lcosto lista;
   typedef lista costos;
   struct tcset{
      char set[Max],aux;
   };
   typedef struct tcset tset;


   unsigned int seed[32];
   int mejor_total[32];
   char perm_total[32][Max];
   int iter_total[32];
   int tproc[32];
   int max_iter,n,optimo;
   int Np;
   int retval;
   unsigned char distancia[Max][Max],flujo[Max][Max];  //matrices de datos


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

    unsigned long mynrand48(unsigned short state[3],unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add)
   {
   
      next(state,mul0,mul1,mul2,add);
      return( ((unsigned long)state[2]) * 0x8000
         + ( ((unsigned long)state[1]) >> 1 )
         );
   
   }


    unsigned long mylrand48(unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in_s[])
   {
      return(mynrand48(in_s,mul0,mul1,mul2,add));
   }



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

// inicializa un conjunto

    void initconta(tset *conjunto)
   {
      conjunto->aux=1;
   }

/*intercambia dos numeros*/
    void cambia(char *a,char *b)
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
// inserta un elemento en la lista 
    void inserta(tlista lista[], unsigned char matriz[][Max])
   {
      int i,j,k=1;
      for(i=1;i<=n;i++)
         for (j=i+1;j<=n;j++)
         {
            lista[k].valor=matriz[i][j];
            lista[k].i=i;
            lista[k].j=j;
            k++;
         }
   }
// intercambia los elementos i,j de la lista 
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
// ordena los valores de la lista tomando en cuenta los indices i,j de su posicion en la matriz
    void ordenaij(tlista a[],int max)
   {
      int j,i;
      for(j=1;j<=max;j++)
         for (i=1;i<=max;i++)
            if (i!=n)
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
//ordena la lista de distancias de forma ascendente
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
// ordena la lista de flujoes de forma descendente
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

// ordena los costos de forma ascendente
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
/* multiplica los flujos grandes y las distancias pequeñas*/
    void mul(costos listacosto[],const int parametro1,tlista lista1[],tlista lista2[])
   {
      int i=1;
      while(i<=parametro1)
      {
         listacosto[i].valor=lista1[i].valor*lista2[i].valor;
         listacosto[i].i=lista2[i].i;
         listacosto[i].j=lista2[i].j;
         listacosto[i].k=lista1[i].i;
         listacosto[i].l=lista1[i].j;
         i++;
      }
      sortcostos(listacosto,1,parametro1);
   
   
   }

// Asigna de forma aleatoria las primeras 2 asignaciones a la permutacion gamma
    void randomiza(const costos listacosto[],parejas gamma[],tset *city,tset *facility,tset *coml,tset *comj,const int parametro2,
    			unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in[])
   {
      int i;
      do
      {
         i=(int)mylrand48(mul0,mul1,mul2,add,in)%(parametro2)+1;
      }
      while (i<=1&&i>=parametro2+1);
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

/* Metodo de ordenacion quicksort*/
    void sortc_gamma(indices a[],int l,int r)
   {
      int aux1,aux2,i,j,x,y;
      i=l;
      j=r;
      x=a[(l+r)/2].valor;
      do {
         while (a[i].valor < x) i++;
         while (x < a[j].valor) j--;
         if (i<=j)
         {
            y=a[i].valor;
            aux1=a[i].j;
            aux2=a[i].i;
            a[i].valor=a[j].valor;
            a[i].j=a[j].j;
            a[i].i=a[j].i;
            a[j].valor=y;
            a[j].j=aux1;
            a[j].i=aux2;
            i++;
            j--;
         }
      }while (i<=j);
      if (l<j) sortc_gamma(a,l,j);
      if (i<r) sortc_gamma(a,i,r);
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
      return t /2;
   }
/* Funcion que devuelve el numero de permutaciones que tienen el mismo costo*/
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


/* funcion que intercambia dos registros usada por los diferentes metodos
de ordenacion*/
    void xchange(indices *a,indices *b)
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


/* ordena de menor a mayor el arreglo c_gamma*/
    void ordena(indices *a[maxlista])
   {
      int j,i;
      for (j=1;j<=n;j++)
         for (i=1;i<=n;i++)
            if (i!=n)
               if (a[i]->valor==a[i+1]->valor)
                  if (a[i]->i==a[i+1]->i)
                     if (a[i]->j>a[i+1]->j)
                        xchange(a[i],a[i+1]);
                     else
                        if (a[i]->i>a[i+1]->i)
                           xchange(a[i],a[i+1]);
   }

/* Construye la permutacion 2 etapa de GRASP*/
    void stage2(parejas gamma[Max], tset *city,tset *facility,tset *coml,tset *comj,
    		unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in[])
   
   {
      int num,aux,suma;
      int i,k,j,rep;
      int m,s;
      indices c_gamma[10000];
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
                     suma=suma+((flujo[i][gamma[j].j])*(distancia[k][gamma[j].l]));
                  aux++;
                  c_gamma[aux].valor=suma;
                  c_gamma[aux].i=i;
                  c_gamma[aux].j=k;
                  m++;
                  suma=0;
               }
         }
         sortc_gamma(c_gamma,1,m);
         s=(int)mylrand48(mul0,mul1,mul2,add,in)%(desempate(c_gamma[1].valor,c_gamma))+1;
         num++;
         addset(facility,c_gamma[s].i);
         outset(comj,c_gamma[s].i);
         addset(city,c_gamma[s].j);
         outset(coml,c_gamma[s].j);
         gamma[num].j=c_gamma[s].i;
         gamma[num].l=c_gamma[s].j;
      }
   }
// calcula los vecinos de gamma con la estructura vecinal lambda
    int lambda(parejas gamma[], const int costo,
    			 unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in[])
   
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
      for (i=1;i<=n;i++)
      {
         addset(&posibles,i);
      }
      do
      {
         j=1;
         k=j+1;
         menor=false;
         do
         {
            if (k==n)
            {
               j++;
               k=j+1;
            }
            if ((inset(&posibles,j)&&!inset(&mov,j)&&!inset(&mov,k))
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
                     iguales[conta].gamma[l].l=gamma[l].l;
                  }
                  conta++;
               }
               cambia(&gamma[k].l,&gamma[j].l);
            }
            k++;
         }while(j!=n-1);
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
         if (j==n-1&&k==n&&menor==false)
            break;
         else
            if (menor!=false)
            {
               for(i=1;i<=n;i++)
               {
                  gamma[i].j=iguales[l].gamma[i].j;
                  gamma[i].l=iguales[l].gamma[i].l;
               }
            }
      }while(1&&error==false);
      if (error!=true)
      {
         for (i=1;i<=n;i++)
         {
            gamma[i].j=iguales[l].gamma[i].j;
            gamma[i].l=iguales[l].gamma[i].l;
         }
         return iguales[l].costo;
      }
      else {
         return aux;
      }
   }
// Calcula los vecinos adyacentes de gamma
    int adyacente(parejas gamma[], const int costo,
    			 unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in[]){
      int aux,total;
      int i, j;
      for (i=1;i<n;i++){
         cambia(&gamma[i].l,&gamma[i+1].l);
         total=costo_total(gamma);
         printgamma(gamma,total);
         if (total<aux) {
            aux=total;
            return total;
         }
         cambia(&gamma[i+1].l,&gamma[i].l);
      
      
      }
      return aux;
   }


/* Calcula los vecinos de gamma con 2 intercambio*/
    int two_exchange(parejas gamma[], const int costo,
    			 unsigned short mul0,unsigned short mul1,
    			unsigned short mul2,unsigned short add,unsigned short in[])
   {
      int i,j,i0,j0,moneda;
      int aux,total;
   
   
      aux = costo;
   
      for (i=1;i<n;i++)
         for (j=i+1;j<=n;j++){
            cambia(&gamma[i].l,&gamma[j].l);
            total=costo_total(gamma);
            if (total<aux) {
               aux=total;
               return total;
            }
            cambia(&gamma[j].l,&gamma[i].l);
         }
        
      return aux;
   }

/*inicializa los conjuntos coml y comj*/
    void init_conjunto(tset *city,tset *facility,tset *coml,tset *comj)
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

// lee el archivo con los datos de las matrices de flujos y de distancias
    void lee_archivo(const char archivo[])
   {
      int i,j;
      float dato;
      FILE *arch;
      if ((arch=fopen(archivo,"r"))==NULL)
      {
         printf("Error de archivo");
         exit(0);
      }
      for (i=1;i<=n;i++){
         for (j=1;j<=n;j++){
            fscanf(arch,"%d",&distancia[i][j]);
         }
      }
      for (i=1;i<=n;i++){
         for (j=1;j<=n;j++){
            fscanf(arch,"%d",&flujo[i][j]);
         }
      }
   }


/* ordena el arreglo de menor a mayor tomando el campo l*/
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
               gamma[j].l=gamma[j+1].l;
               gamma[j+1].j=aux.j;
               gamma[j+1].l=aux.l;
            }
   }




/* inicializa la permutacion poniendo 0*/
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
         destino[i].j=fuente[i].j;
         destino[i].l=fuente[i].l;
      }
   }
   
    void grasp()
   {
      static unsigned short internal_state[3] = {1, 0, 0};
      static unsigned short multiplier0 = RAND48_MULT0;
      static unsigned short multiplier1 = RAND48_MULT1;
      static unsigned short multiplier2 = RAND48_MULT2;
      static unsigned short additiver = RAND48_ADD;
   // Variables para tomar el tiempo de ejecucion
      //struct rusage ru;
      float t_inicial,t_final;
      float timep;
      time_t t;
   
      struct tcset city,facility,coml,comj,mov,posibles; // variables para conjuntos
      parejas gamma[Max];  //permutacion
      float beta,alpha; 	// parametros para restringuir la cantidad de elementos de las listas de flujos y de distancias
      tlista lflujo[maxlista],ldist[maxlista];  // lista de flujos y de distancias
      int total,res,parametro1,parametro2,band,menor;
      char archivo[8],archivo1[12]="sal";
      costos listacosto[maxlista]; // lista de costos
      int best,iter_local; //
      unsigned long  iter_grasp,iter_mejor;
      parejas mejorp[Max];
      int mejorcosto = 32000;
      int i,j,my_ide = 0;
      FILE *outfile;
      
   	// Se toma el tiempo inicial
      getrusage(RUSAGE_SELF,&ru);
      t_inicial = (float)ru.ru_utime.tv_sec +
         	 (float)(ru.ru_utime.tv_usec)/1000000;
      mysrand48(seed[my_ide],&multiplier0,&multiplier1,&multiplier2,&additiver,internal_state);
   // inicializa los contadores de iteracion
      iter_grasp=iter_mejor=0;
      iter_local=0;
      
   	// establece los parametros alpha y beta
      res=(n*(n-1))/2;
      beta=.1;
      parametro1=(int)(res*beta)+1;
      alpha=.5;
      parametro2=(int)(parametro1*alpha)+1;
      
   	// inserta los elementos de la matriz de flujos en la lista de flujos
      inserta(lflujo,flujo);
      // ordena la lista de flujos
      sortlflujo(lflujo,1,res);
      ordenaij(lflujo,res);
      
   	// inserta los elementos de la matriz de distancias en la lista de distancias
      inserta(ldist,distancia);
       // ordena la lista de distancias
      sortldist(ldist,1,res);
      ordenaij(ldist,res);
      // multiplica la lista de flujos por la de distancias 
      mul(listacosto,parametro1,ldist,lflujo);
      
   	// inicializa la lista de movimientos posibles
      initconta(&city);
      initconta(&facility);
      initconta(&coml);
      initconta(&comj);
   // repite hasta alcanzar el maximo numero de iteraciones o alcanzar o mejorar el valor optimo o mejor valor conocido
      do
      {
      // inicializa la permutacion gamma
         limpiagamma(gamma);
         // inicializa los conjuntos de movimientos con valores validos
         init_conjunto(&city,&facility,&coml,&comj);
         // asigna los primeros 2 elementos de la permutacion
         randomiza(listacosto,gamma,&city,&facility,&coml,&comj,parametro2,multiplier0,multiplier1,multiplier2,additiver,internal_state);
         // calcula y asigna los n - 2 elementos restantes a la permutacion
         stage2(gamma, &city,&facility,&coml,&comj,multiplier0,multiplier1,multiplier2,additiver,internal_state);
         // ordena la permutacion
         buble(gamma);
      	// calcula su costo inicial
         best=costo_total(gamma);
         iter_local=0;
         menor =true;
         // repite hasta que ya no se encuentre una mejor permutacion
         while(menor==true){
            menor=false;
            // repite la estructura vecinal  hasta que no se encuentre una mejor permutacion
            do  
            {
               total=adyacente(gamma, best,multiplier0,multiplier1,multiplier2,additiver,internal_state);
               if (total<best)
               {
                  band=true;
                  best=total;
                  menor=true;
                  copiap(gamma,mejorp);
                  mejorcosto=best;
               
               }
               else
                  band=false;
            }while (band!=false);
            band=false;
            // repite la estructura vecinal  hasta que no se encuentre una mejor permutacion
         
            do  
            {
               total=two_exchange(gamma, best,multiplier0,multiplier1,multiplier2,additiver,internal_state);
               printgamma(gamma,total);
               if (total<best)
               {
                  band=true;
                  best=total;
                  menor=true;
                  copiap(gamma,mejorp);
                  mejorcosto=best;
               
                 
               }
               else
                  band=false;
            }while (band!=false);
            band=false;
           // repite la estructura vecinal  hasta que no se encuentre una mejor permutacion
         
            do  
            {
               total=lambda(gamma, best,multiplier0,multiplier1,multiplier2,additiver,internal_state);
               printgamma(gamma,total);
               if (total<best)
               {
                  band=true;
                  best=total;
                  menor=true;
                  copiap(gamma,mejorp);
                  mejorcosto=best;
               
               }
               else
                  band=false;
            }while (band!=false);
            band=false;
            iter_local++;
         }
      	// si se encontro un mejor valor se guarda
         if (best<mejorcosto){
            copiap(gamma,mejorp);
            mejorcosto=best;
            iter_mejor = iter_grasp;
         }
         iter_grasp++;
      }while(((iter_grasp-iter_mejor)<max_iter)&&(optimo<mejorcosto));
      
   	// se toma el tiempo final
      getrusage(RUSAGE_SELF,&ru);
      t_final = (float)ru.ru_utime.tv_sec +
         	 (float)(ru.ru_utime.tv_usec)/1000000;
         	 
   	// se calcula el tiempo de procesamiento			 
      timep = (t_final- t_inicial);
   
      mejor_total[my_ide] = mejorcosto;
      iter_total[my_ide] = iter_grasp;
      for (i=1;i<=n;i++){
         perm_total[my_ide][i] =  mejorp[i].l;
      }
      tproc[my_ide]=timep;
   
   // se almacena el mejor valor econtrado en el archivo de salida
      sprintf(archivo,"%d",my_ide);
      strcat(archivo1,archivo);
      outfile = fopen (archivo1,"w");
      fprintf (outfile,"my_id = %7d\n",my_ide);
      for (j=1;j<=n;j++){
      
         fprintf (outfile,"%3d",perm_total[my_ide][j]);
      }
      fprintf (outfile,"\n");
      fprintf (outfile,"%7d",mejor_total[my_ide]);
      fprintf (outfile,"%7d",iter_total[my_ide]);
      fprintf (outfile,"%12.4f",tproc[my_ide]);
      fclose (outfile);
   
   
   }

    void vacia(void){
   }

    void main(argc, argv)
       int argc;
       char *argv[];
   {
      FILE *outfile;
      //struct rusage ru;
      float t_inicial,t_final;
      float timep;
      time_t t;
      int i,j;
   
   
      if (argc != 7){
         printf("Forma de uso:\n");
         printf("Archivo de datos\n");
         printf("numero del problema\n");
         printf("numero maximo de iteraciones\n");
         printf("optimo\n");
        // printf("numero de procesadores\n");
         printf("archivo de salida\n");
      
         exit (1213);
      }
   
      n = atoi(argv[2]); // instancia del problema
      max_iter = atoi(argv[3]);  // numero maximo de iteraciones
      optimo = atoi(argv[4]);   // optimo o mejor valor conocido como meta 
     //Np = atoi(argv[5]);
   
      outfile = fopen (argv[5],"w");
   //m_set_procs(Np);
      fprintf (outfile,"Np = %d\n",Np);
   
      for (i=0;i<Np;i++)
         seed[i] = time(&t)*(i+1);
   
   // lectura del archivo de datos
      lee_archivo(argv[1]);
   
      getrusage(RUSAGE_SELF,&ru);
      t_inicial = (float)ru.ru_utime.tv_sec +
         	 (float)(ru.ru_utime.tv_usec)/1000000;
   
   //ejecucion del grasp
      grasp();
   	   
      getrusage(RUSAGE_SELF,&ru);
      t_final = (float)ru.ru_utime.tv_sec +
         	 (float)(ru.ru_utime.tv_usec)/1000000;
      timep = (t_final- t_inicial);
   
      fprintf (outfile,"Np = %d\n",Np);
      fprintf (outfile,"tiempo %12.4f\n",timep);
   
      for (i=0;i<Np;i++){
         for (j=1;j<=n;j++)
            fprintf (outfile,"%3d",perm_total[i][j]);
         fprintf (outfile,"%7d",mejor_total[i]);
         fprintf (outfile,"%7d",iter_total[i]);
         fprintf (outfile,"%12.4f",tproc[i]);
         fprintf (outfile,"\n");
      }
      fclose (outfile);
   
   }
