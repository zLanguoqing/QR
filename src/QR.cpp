#include "QR.h"
/* 
 一般实矩阵的QR分解
 采用householder进行QR迭代
 参数说明
 a 存放实矩阵的数组，返回时其右上三角存放的是QR分解的R
 m 矩阵行数  n 矩阵列数
 q 返回时正交矩阵Q
*/
int QR(double a[], int m, int n, double q[])
{
    int i, j, k, l, nn, p, jj;
    double u, alpha, w, t;
    if (m < n)
    {
        printf("fail\n");
        return 0;
    }
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= m - 1; j++)
        {
            l = i * m + j;
            q[l] = 0.0;
            if (i == j)
                q[l] = 1.0;
        }
    }
    nn = n;
    if (m == n)
    {
        nn = m - 1;
    }
    for (k = 0; k <= nn - 1; k++)
    {
        u = 0.0;
        l = k * n + k;
        for (i = k; i <= m - 1; i++)
        {
            w = fabs(a[i * n + k]);
            if (w > u)
                u = w;
        }
        alpha = 0.0;
        for (i = k; i <= m - 1; i++)
        {
            t = a[i * n + k] / u;
            alpha = alpha + t * t;
        }
        if (a[l] > 0.0)
            u = -u;
        alpha = u * sqrt(alpha);
        if (fabs(alpha) + 1.0 == 1.0)
        {
            printf("fail\n");
            return 0;
        }
        u = sqrt(2.0 * alpha * (alpha - a[l]));
        if (u > MIN)
        {
            a[l] = (a[l] - alpha) / u;
            for (i = k + 1; i <= m - 1; i++)
            {
                p = i * n + k;
                a[p] = a[p] / u;
            }
            for (j = 0; j <= m - 1; j++)
            {
                t = 0.0;
                for (jj = k; jj <= m - 1; jj++)
                {
                    t = t + a[jj * n + k] * q[jj * m + j];
                }
                for (i = k; i <= m - 1; i++)
                {
                    p = i * m + j;
                    q[p] = q[p] - 2.0 * t * a[i * n + k];
                }
            }
            for (j = k + 1; j <= n - 1; j++)
            {
                t = 0.0;
                for (jj = k; jj <= m - 1; jj++)
                {
                    t = t + a[jj * n + k] * a[jj * n + j];
                }
                for (i = k; i <= m - 1; i++)
                {
                    p = i * n + j;
                    a[p] = a[p] - 2.0 * t * a[i * n + k];
                }
            }
            a[l] = alpha;
            for (i = k + 1; i <= m - 1; i++)
            {
                a[i * n + k] = 0.0;
            }
        }
    }
    for (i = 0; i <= m - 2; i++)
    {
        for (j = i + 1; j <= m - 1; j++)
        {
            p = i * m + j;
            l = j * m + i;
            t = q[p];
            q[p] = q[l];
            q[l] = t;
        }
    }
    return 1;
}

/*
将一般实矩阵约化为Hessenberg矩阵
 参数说明
 a 存放一般实矩阵，返回为上H矩阵
 n 矩阵的阶数
*/
void Hessenberg(double a[], int n)
{
    int i, j, k, u, v;
    double d, t;
    for (k = 1; k <= n - 2; k++)
    {
        d = 0.0;
        for (j = k; j <= n - 1; j++)
        {
            u = j * n + k - 1;
            t = a[u];
            if (fabs(t) > fabs(d))
            {
                d = t;
                i = j;
            }
        }
        if (fabs(d) > MIN)
        {
            if (i != k)
            {
                for (j = k - 1; j <= n - 1; j++)
                {
                    u = i * n + j;
                    v = k * n + j;
                    t = a[u];
                    a[u] = a[v];
                    a[v] = t;
                }
                for (j = 0; j <= n - 1; j++)
                {
                    u = j * n + i;
                    v = j * n + k;
                    t = a[u];
                    a[u] = a[v];
                    a[v] = t;
                }
            }
            for (i = k + 1; i <= n - 1; i++)
            {
                u = i * n + k - 1;
                t = a[u] / d;
                a[u] = 0.0;
                for (j = k; j <= n - 1; j++)
                {
                    v = i * n + j;
                    a[v] = a[v] - t * a[k * n + j];
                }
                for (j = 0; j <= n - 1; j++)
                {
                    v = j * n + k;
                    a[v] = a[v] + t * a[j * n + i];
                }
            }
        }
    }
    return;
}
/*
a 存放上H矩阵
n 上H矩阵的阶数
u 返回n个特征值的实部
v 返回n个特征值的虚部
esp 控制精度要求
Iteration 迭代次数
*/
int HessenbergQR(double a[], int n, double u[], double v[], double eps, int Iteration)
{
    int m, it, i, j, k, l, ii, jj, kk, ll;
    double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
    it = 0;
    m = n;
    while (m != 0)
    {
        l = m - 1;
        while ((l > 0) && (fabs(a[l * n + l - 1]) > eps *
                                                        (fabs(a[(l - 1) * n + l - 1]) + fabs(a[l * n + 1]))))
        {

            l = l - 1;
        }
        ii = (m - 1) * n + m - 1;
        jj = (m - 1) * n + m - 2;
        kk = (m - 2) * n + m - 1;
        ll = (m - 2) * n + m - 2;
        if (l == m - 1)
        {
            u[m - 1] = a[(m - 1) * n + m - 1];
            v[m - 1] = 0.0;
            m = m - 1;
            it = 0;
        }
        else
        {
            if (l == m - 2)
            {
                b = -(a[ii] + a[ll]);
                c = a[ii] * a[ll] - a[jj] * a[kk];
                w = b * b - 4.0 * c;
                y = sqrt(fabs(w));
                if (w > MIN)
                {
                    xy = 1.0;
                    if (b < MIN)
                        xy = -1.0;
                    u[m - 1] = (-b - xy * y) / 2.0;
                    u[m - 2] = c / u[m - 1];
                    v[m - 1] = 0.0;
                    v[m - 2] = 0.0;
                }
                else
                {
                    u[m - 1] = -b / 2.0;
                    u[m - 2] = u[m - 1];
                    v[m - 1] = y / 2.0;
                    v[m - 2] = -v[m - 1];
                }
                m = m - 2;
                it = 0;
            }
            else
            {
                if (it >= Iteration)
                {
                    printf("fail \n");
                    return -1;
                }
                it = it + 1;
                for (j = l + 2; j <= m - 1; j++)
                    a[j * n + j - 2] = 0.0;
                for (j = l + 3; j <= m - 1; j++)
                    a[j * n + j - 3] = 0.0;
                for (k = l; k <= m - 2; k++)
                {
                    if (k != l)
                    {
                        p = a[k * n + k - 1];
                        q = a[(k + 1) * n + k - 1];
                        r = 0.0;
                        if (k != m - 2)
                            r = a[(k + 2) * n + k - 1];
                    }
                    else
                    {
                        x = a[ii] + a[ll];
                        y = a[ll] * a[ii] - a[kk] * a[jj];
                        ii = l * n + l;
                        jj = l * n + l + 1;
                        kk = (l + 1) * n + l;
                        ll = (l + 1) * n + l + 1;
                        p = a[ii] * (a[ii] - x) + a[jj] * a[kk] + y;
                        q = a[kk] * (a[ii] + a[ll] - x);
                        r = a[kk] * a[(l + 2) * n + l + 1];
                    }
                    if ((fabs(p) + fabs(q) + fabs(r)) > MIN)
                    {
                        xy = 1.0;
                        if (p < MIN)
                            xy = -1.0;
                        s = xy * sqrt(p * p + q * q + r * r);
                        if (k != l)
                            a[k * n + k - 1] = -s;
                        e = -q / s;
                        f = -r / s;
                        x = -p / s;
                        y = -x - f * r / (p + s);
                        g = e * r / (p + s);
                        z = -x - e * q / (p + s);
                        for (j = k; j <= m - 1; j++)
                        {
                            ii = k * n + j;
                            jj = (k + 1) * n + j;
                            p = x * a[ii] + e * a[jj];
                            q = e * a[ii] + y * a[jj];
                            r = f * a[ii] + g * a[jj];
                            if (k != m - 2)
                            {
                                kk = (k + 2) * n + j;
                                p = p + f * a[kk];
                                q = q + g * a[kk];
                                r = r + z * a[kk];
                                a[kk] = r;
                            }
                            a[jj] = q;
                            a[ii] = p;
                        }
                        j = k + 3;
                        if (j > m - 1)
                            j = m - 1;
                        for (i = l; i <= j; i++)
                        {
                            ii = i * n + k;
                            jj = i * n + k + 1;
                            p = x * a[ii] + e * a[jj];
                            q = e * a[ii] + y * a[jj];
                            r = f * a[ii] + g * a[jj];
                            if (k != m - 2)
                            {
                                kk = i * n + k + 2;
                                p = p + f * a[kk];
                                q = q + g * a[kk];
                                r = r + z * a[kk];
                                a[kk] = r;
                            }
                            a[jj] = q;
                            a[ii] = p;
                        }
                    }
                }
            }
        }
    }
    return 1;
}
/*
矩阵相乘
 a a矩阵
 b b矩阵
 m a矩阵与两个矩阵乘积c的行数
 n a矩阵的列数与b矩阵的列数
 k 矩阵b与c的列数
 c 返回 c = ab
*/
void Mul(double a[], double b[], int m, int n, int k, double c[])
{
    int i, j, l, u;
    for (i = 0; i <= m - 1; i++)
        for (j = 0; j <= k - 1; j++)
        {
            u = i * k + j;
            c[u] = 0.0;
            for (l = 0; l <= n - 1; l++)
                c[u] = c[u] + a[i * n + l] * b[l * k + j];
        }
}
/*
矩阵求逆
采用G-L消去法
a 输入矩阵
n 矩阵的阶数 
*/
int Inv(double a[], int n)
{
    int *is, *js, i, j, k, l, u, v;
    double d, p;
    is = (int *)malloc(n * sizeof(int));
    js = (int *)malloc(n * sizeof(int));
    for (k = 0; k <= n - 1; k++)
    {
        d = 0.0;
        for (i = k; i <= n - 1; i++)
            for (j = k; j <= n - 1; j++)
            {
                l = i * n + j;
                p = fabs(a[l]);

                if (p > d)
                {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }

        if (d +1.0 == 1.0)
        {
            free(is);
            free(js);
            printf("err * * not inv\n");
            return 0;
        }
        if (is[k] != k)
            for (j = 0; j <= n - 1; j++)
            {
                u = k * n + j;
                v = is[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        if (js[k] != k)
            for (i = 0; i <= n - 1; i++)
            {
                u = i * n + k;
                v = i * n + js[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        l = k * n + k;
        a[l] = 1.0 / a[l];
        for (j = 0; j <= n - 1; j++)
            if (j != k)
            {
                u = k * n + j;
                a[u] = a[u] * a[l];
            }
        for (i = 0; i <= n - 1; i++)
            if (i != k)
                for (j = 0; j <= n - 1; j++)
                    if (j != k)
                    {
                        u = i * n + j;
                        a[u] -=  a[i * n + k] * a[k * n + j];
                    }
        for (i = 0; i <= n - 1; i++)
            if (i != k)
            {
                u = i * n + k;
                a[u] = -a[u] * a[l];
            }
    }
    for (k = n - 1; k >= 0; k--)
    {
        if (js[k] != k)
            for (j = 0; j <= n - 1; j++)
            {
                u = k * n + j;
                v = js[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        if (is[k] != k)
            for (i = 0; i <= n - 1; i++)
            {
                u = i * n + k;
                v = i * n + is[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
    }
}
