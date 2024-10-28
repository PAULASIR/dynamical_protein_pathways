/* VDP oscillator as a model of PPI */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>


#define MAX_SIZE 121
#define BUFFER_SIZE 601


int i,j,k,nn=100,n=2;
double mu;
int rc, rr;
double y[MAX_SIZE][BUFFER_SIZE];
double z[250000], zz[250000];

double arr[101][101], arr1[101][101], arr2[101][101], a[101][101];
double sm, diff, eps, eps1, eps2;
double complex  phi1, phi2;
int pc;
double p;
double syn_err, sm, avg;
double mean[MAX_SIZE];

double sum;

double summ[MAX_SIZE], *num;



void RK4(int,int,double,double,double[MAX_SIZE][BUFFER_SIZE],
                   void (*DGL)(double,double[MAX_SIZE][BUFFER_SIZE],double[MAX_SIZE][BUFFER_SIZE]));
void DGL(double, double[MAX_SIZE][BUFFER_SIZE],double[MAX_SIZE][BUFFER_SIZE]);


/****************************function to generate random coupling matrices*********/
void adj_mat(int row, int col, int arr[row][col], int arr1[row][col], int arr2[row][col])
{
    
  //  srand(time(0));
    
    for (int i = 0; i < row; i++)
      {
        for (int j = 0;  j < col; j++) {
           arr2[i][j] = rand()%2;
           if (i <  pc && j< pc)
                {arr[i][j] = 1;}
             else
                {arr[i][j] = 0;}
                
           
            rc = rand()%nn;
            rr= rand()%nn;
              
        int temp =  arr[i][j];
         arr[i][j] =  arr[rr][rc];
         arr[rr][rc] = temp;
              
  
            if (arr[i][j] == 1)
               { arr1[i][j]  = 0; }
                 else
               { arr1[i][j]  = 1; } 
           }   
           //  printf(" \n");
    }
  }
/*************************************/  
   double *getarray(double *summ)  
    {  
        for(int i=1;i<=nn;i++)  
          {  
             &summ[i] ; 
          }  
        return summ;  
    } 
/**********************************/
void main()
 {
//nn= no 0f oscillators, n=dimension of the model//
   double t,h;
  
  // srand(time(0));

FILE *fp1,*fp2,*fp3,*fp4;
fp1=fopen("syn_rate_eps.dat","w");


for(pc = 0; pc <= 100; pc+= 2) {


//printf("%d \n", pc);
 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0+2.0;
 }




eps1 = 0.93; eps2 = 0.93; eps = 0.93;
// PARAMETERS
mu = 1.0;
 
//***time step***//
h=0.01; t=0.0;



for(k=1;k<=3000;k++)
   {               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);   
    	
   for(j=1;j<=nn;j++)
      {  
/*    leaving transients*/               
     if(k>=10)
	{                                                           
           z[k] = y[j][1]; 
       
                            
        phi1=cexp(y[1][1]*I);   phi2=cexp(y[j-1][1]*I);    
        phi1=carg(phi1);        phi2=carg(phi2);    
        
      // printf("%f \n", creal(phi1) - creal(phi2) );
        diff = fabs(creal(phi1) - creal(phi2));  
        
          if (diff  < 0.00001 )
             summ[j] =  k;
             
          }
     }    
}
 
    num = getarray(summ);
  
     sum=0.0;
     for(int i=1; i<= nn; i++)  
        {  
/*           if(num[i] != num[i])*/
/*              { num[i] = 0;}*/
/*            else */
/*              { num[i] = num[i]; }*/
/*               */
                 sum += num[i];  

                 // fprintf(fp1,"%f   %f\n", p,num[i]);

        }  
           
         
         sum = sum/(float) (3000*nn);
           
           printf("%f \n", sum);
           
           double p = pc/(float) nn;
         
               fprintf(fp1,"%f   %f\n", p, sum );
              
              
               //fprintf(fp1,"%f   %f\n", eps,num[i]);


/*    for(i =1; i<=nn; i++){*/
/*                mean[i] = sum/ nn;*/
/*                avg = (fabs((mean[i])) + fabs((mean[i+1])))/ 2.0 ;}*/

    
          // fprintf(fp1,"%f  %f  %f\n", p, avg, avg/ (2*M_PI));
 }
/******************************************************/
    
  printf("process over!!\n");
}
//************************RK4 SUBROUTINE*********************************//
void RK4(int nn,int n,double h,double t,double y[MAX_SIZE][BUFFER_SIZE],
	   void (*DGL)(double,double[MAX_SIZE][BUFFER_SIZE],double[MAX_SIZE][BUFFER_SIZE]))
{
	   
	   double k1[MAX_SIZE][BUFFER_SIZE],k2[MAX_SIZE][BUFFER_SIZE],k3[MAX_SIZE][BUFFER_SIZE],k4[MAX_SIZE][BUFFER_SIZE];
	   double yaux[MAX_SIZE][BUFFER_SIZE];

	   DGL(t,y,k1);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k1[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k2);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k2[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k3);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k3[j][i];
	   }
	   
	   DGL(t+h,yaux,k4);
	   for(j=1;j<=nn;j++)
	   {
             for(i=1;i<=n;i++)
	      y[j][i]=y[j][i]+h*((k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.0);
	   }
}
//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[MAX_SIZE][BUFFER_SIZE],double F[MAX_SIZE][BUFFER_SIZE])
 {
     int row = nn;
     int col = nn;
     int arr[row][col];
     int arr1[row][col];
     int arr2[row][col];
     
      adj_mat(row, col, arr, arr1, arr2);
      
     /*******************************************************/

   for(j=1;j<=nn;j++)
     { 
       for(i=1;i<=nn;i++)
          {
         //  printf("%d   %d  %d\n", arr[i][j], arr1[i][j], arr2[i][j]);
       F[j][1]= y[j][2] + eps*arr[i][j]*(y[i][1]-y[j][1]) + eps1*arr1[i][j]*(y[i][1]-y[j][1]) + eps2*arr2[i][j]*(y[i][1]-y[j][1]); 
       F[j][2]= mu*y[j][2]*(1-y[j][1]*y[j][1]) - y[j][1];      
    }
  }
}
/************calculating mean phase*******

for(j=0; j<=nn; j++){

  double x_bar, y_bar,phi_avg, p_avg;    
  
      y_bar = *(ptr_y+j)/(double)(del_t); 
      x_bar = *(ptr_x+j)/(double)(del_t); 
           

        printf("%d  %f  \n",j, *(ptr_y+j));  

           phi_avg = atan(y_bar/x_bar);

           
           if( (phi_avg > 2.0*M_PI) || (phi_avg<0.0) )
              { p_avg =0.0; } 
              else
                { p_avg = phi_avg;} 
              
         //  fprintf(fp4,"%d  %f\n", j, p_avg); 
        }
/**************************************/




/*----------------max & min of the series--------
       if(y[j][1]>x_max){x_max=y[j][1];} 
       if(y[j][1]<x_min){x_min=y[j][1];}   */


/*****************return y[j][1]*******************
double getsum( double summ[], int nn)
  {
     int i; double sm;
       for(i =1;i<= nn; i++){
           sm = summ[i] ;
              }
                
    return sm;
  }
/*********************************/
