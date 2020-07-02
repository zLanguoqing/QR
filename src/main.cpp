#include "stdio.h"
#include "QR.h"
int main()
{
    int i, j;
    double q[16], a[12] = {1.0, 1.0, -1.0, 2.0, 1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 2.0, 1.0} ;
    i = QR(a, 4, 3, q);
    printf("\n");
	printf("i=%d\n", i);
	printf("\nMAT Q Is:\n");
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
			printf("%e ", q[i * 4 + j]);
		printf("\n");
	}
    printf("\nMAT R Is:\n");
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 3; j++)
			printf("%e ", a[i * 3 + j]);
		printf("\n");
    }
    return 0;
}