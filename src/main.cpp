#include "stdio.h"
#include "QR.h"
int main()
{
    int i, j;
    double q[16], a[12] = {1.0, 1.0, -1.0, 
                            2.0, 1.0, 0.0, 
                            1.0, -1.0, 0.0,
                            -1.0, 2.0, 1.0} ;
    i = QR(a, 4, 3, q);
    printf("\n");
	printf("i=%d\n", i);
	printf("\nMAT Q Is:\n");
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
			printf("%13.6e ", q[i * 4 + j]);
		printf("\n");
	}
    printf("\nMAT R Is:\n");
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 3; j++)
			printf("%13.6e ", a[i * 3 + j]);
		printf("\n");
    }
    int ii,jj;
    double A[25] = {1.0, 6.0, -3.0, -1.0, 7.0,
                    8.0, -15.0, 18.0, 5.0, 4.0,
                    -2.0, 11.0, 9.0, 15.0, 20.0,
                    -13.0, 2.0, 21.0, 30.0, -6.0,
                    17.0, 22.0, -5.0, 3.0, 6.0};
    Hessenberg(A,5);
    printf("MAT A IS:\n");
    for(ii = 0; ii <= 4;ii++)
    {
        for(jj = 0;jj <= 4;jj++)
            printf("%13.6e ",A[ii*5+jj]);
        printf("\n");
    }
    printf("\n");
    double u[5],v[5];
    double eps = 0.000001;
    i = HessenbergQR(A,5,u,v,eps,60);
    if(i > 0)
    {
        for(i = 0; i <= 4; i++)
            printf("%13.6e+J %13.6e \n",u[i],v[i]);
    }
    printf("\n");
    return 0;
}