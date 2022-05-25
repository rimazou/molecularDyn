#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#define NUM_THREADS 8
#define n 100


typedef struct thread_args
{
	int my_first;
	int my_last;
        double seed;
} thread_args;
typedef struct thread_args thread_args;
double r[200];
double a,b;

clock_t start, end;
double temps;
void *r8mat_uniform_ab(void *input)

{
  thread_args *args = (thread_args *)input;
	int my_first = args->my_first;
	int my_last = args->my_last;
    double seed = args->seed;
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }  

  for ( j = 0; j < n; j++ )
  {   
    for (i = my_first; i < my_last; ++i)
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;
       
      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
    
      r[i+j*my_last] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }


}
int main (void){
 
 pthread_t threads[NUM_THREADS];
 int rs;
 int i;
 
   	start = clock();

 	for (i = 0; i < NUM_THREADS; i++)
	{
		thread_args *args = malloc(sizeof(thread_args));
		args->my_first = i * (n / NUM_THREADS);
		args->my_last = (i + 1) * (n / NUM_THREADS) - 1;
        args->seed = 5;
		rs = pthread_create(&threads[i], NULL, r8mat_uniform_ab, (void*)&args);
       
		if (rs)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rs);
			exit(-1);
		}
     
	}
	for (i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
       
	}
 
  
   
  	end = clock();
 
  	temps = ((double)(end - start)) / CLOCKS_PER_SEC;
  	printf ("Calcul parallel %f secondes\n", temps);
   
  	
  	
}
