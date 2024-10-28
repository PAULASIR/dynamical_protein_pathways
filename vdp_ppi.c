/* VDP oscillator as a model of PPI */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>


#define MAX_SIZE 121
#define BUFFER_SIZE 601


int i,j,k,nn=50,n=2;
double mu;
int rc, rr;
double y[MAX_SIZE][BUFFER_SIZE];
double z[250000], zz[250000];
//double *ptr, num[];
double arr[110][110], arr1[110][110], arr2[110][110], a[110][110];
double sm,diff,mean, eps, eps1, eps2;
double p;
int pc;

double summ[MAX_SIZE];


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
                
            rc = rand()%nn-1;
            rr= rand()%nn-1;
              
        int temp =  arr[i][j];
         arr[i][j] =  arr[rr][rc];
         arr[rr][rc] = temp;
              

        //printf("%d" ,arr[i][j]);
            if (arr[i][j] == 1)
               { arr1[i][j]  = 0; }
                 else
               { arr1[i][j]  = 1; } 
           }   
           //  printf(" \n");
    }
  }
/*****************return y[j][1]*********************/
double getsum( double summ[], int nn)
  {
     int i; 
       for(i =1;i<= nn; i++)
          printf("%f \n", summ[i]);
              
    return 0;
  }
/*********************************/
void main()
 {
//nn= no 0f oscillators, n=dimension of the model//
   double t,h;
  
  // srand(time(0));

FILE *fp1,*fp2,*fp3,*fp4;
fp1=fopen("1.dat","w");
fp2=fopen("2.dat","w");
fp3=fopen("3.dat","w");
fp4=fopen("4.dat","w");


 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0+2.0;
 }


//PARAMETERS
printf("eps1=\n");
scanf("%lf",&eps1);
printf("eps=\n");
scanf("%lf",&eps);
printf("eps2=\n");
scanf("%lf",&eps2);


mu = 2.0;
 
//***time step***//
h=0.08; t=0.0;

pc=20;

double p = pc/(float) nn;

for(k=1;k<=800;k++)
  {               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);   
    	
   for(j=1;j<=nn;j++)
      {  
/*    leaving transients*/               
	if(k>=10)
	{                                                           
          z[k] = y[j][1];      
                   
/* -----------------Printing time series-----------------*/        
       if( k == 500){
              //summ[j] = y[j][1];
           fprintf(fp3,"%d %f \n",j,y[j][1]);} 
   
      fprintf(fp2,"%d %f %f  %f\n",j,t,y[j][1],y[j][2]); 
    
       fprintf(fp1,"%f\t",t);
         for(i=1;i<=nn;i++)
           fprintf(fp1,"  %f\t",y[i][1]);    
           fprintf(fp1,"\n");
/**************************************************/    
          }
     }    
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
  
    /*****************************************************  
     
    meanx=0.0; 

    for(j=1;j<=nn;j++)
     {
       meanx=meanx+y[j][1];
     }

    meanx=meanx/nn; 
     /*******************************************************/

   for(j=1;j<=nn;j++)
     { 
       for(i=1;i<=nn;i++)
          {
   //       printf("%d   \n",  arr[i][j]);
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

