    #include <stdio.h>  
    int *getarray(int *a)  
    {  
        
        printf("Enter the elements in an array : ");  
        for(int i=0;i<5;i++)  
          {  
          //  scanf("%d", &a[i]);
          &a[i] ;
            
          }  
        return a;  
    }  
    
    int main()  
     {  
      int *n;  
      int a[5]; 
      
      for(int i =0; i<5; i++){
         a[i] = i;}
      
       
      n=getarray(a);  
      printf("\n Elements of array are :");  
      for(int i=0; i<5; i++)  
        {  
            printf("%d", n[i]);  
        }  
        return 0;  
     }  
