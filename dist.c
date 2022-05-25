#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#define SIZE 3
#define NUM_THREADS 3
#define CHUNK SIZE/NUM_THREADS

int id[NUM_THREADS];
clock_t start, end;
double temps;

double d, r1[SIZE], r2[SIZE], dr[SIZE];
pthread_t tid[NUM_THREADS];
pthread_mutex_t mutex_d;

void* dist(void* id) {
  size_t i;
  int my_first = *(int*)id * CHUNK;
  int my_last = (*(int*)id + 1) * CHUNK;
  double loc = 0.;
  
  // Calcul
  for (i = my_first; i < my_last; i++)
    {
    dr[i] = r1[i] - r2[i];
    loc = loc + dr[i]*dr[i];
    
  }
  pthread_mutex_lock(&mutex_d);
  d = d + loc;
  pthread_mutex_unlock(&mutex_d);
  return NULL;
}

int main() {
  size_t i;
  
  // Initialisation
  d = 0.;
  for (i = 0; i < SIZE; i++) {
    r1[i] = i * 0.4;
    r2[i] = i * 3.0;
  }
  
 pthread_mutex_init(&mutex_d, NULL);
  start = clock();
  for (i = 0; i < NUM_THREADS; i++) {
    id[i] = i;
    pthread_create(&tid[i], NULL, dist,
                   (void*)&id[i]);
  }
  
  for (i = 0; i < NUM_THREADS; i++)
  pthread_join(tid[i], NULL);
  end = clock();
  pthread_mutex_destroy(&mutex_d);

  printf("d = %g\n", sqrt(d));
    
  temps = ((double)(end - start)) / CLOCKS_PER_SEC; 
  printf ("Calcul parallèle %f secondes\n", temps);
  return 0;
}
