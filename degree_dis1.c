#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>


void swap (int *a, int *b)
 {
    int temp = *a;
    *a = *b;
    *b = temp;
 }s


/*-----------------------------------------------------------*/
int randomize(int aa[], int nn)
  {
    srand ( time(NULL));
    for (int i = nn-1; i > 0; i--){
        int j = rand() % (i+1);
        swap(&aa[i], &aa[j]);
    }
  }
/*--------------------------------------------------------*/

void main()
 {

   int nn=10; 
   int i,j;
   int arr[nn][nn], a[nn][nn], index[nn];
   int deg_sum;
   int n=100;
   int pc = 50;
   double p = pc/(float) n;

   printf("%f\n", p);
  

 int shuffle(int a[][nn], int nn){
 int i,j;
  for (int i= 1; i < nn; i++){      
        int r = rand()%nn+1;
      // int r = rand() % (i+1);
       //int r = (rand() % (nn - i)) + i;
       //  printf("%d\n", r);
     for (int j = 1; j < nn; j++)
       {
         int temp =  a[i][j];
         a[i][j] =  a[r][i];
         a[r][j] = temp;
       }
  }
}
  
   
   srand(time(0)); 
   deg_sum =0.0;
    
    for(int i = 0; i < 10; i++){
       for(int j = 0; j < 10; j++){
             if (j <= pc && i <= pc)
                arr[i][j] = 1;
             else
                arr[i][j] =0;
                         
                        
             randomize(arr,nn);     
             deg_sum += arr[i][j];        
            printf("%d  " ,arr[i][j] );
        }
        
        printf("\n");
       } 
    

        printf("%d \n", deg_sum);

}






