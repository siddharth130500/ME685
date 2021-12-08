#include<bits/stdc++.h>
using namespace std;

double f(double x)
{
    return x*exp(-x);  //The function to be integrated
}

double trapezoidal(int N)
{
    double deltax=(double)(3-0)/(double)(N-1); //Dividing the region into N-1 equally spaced points
    vector<double> x(N);    //An array that includes the value of x at every node
    x[0]=0;
    for(int i=1;i<N;i++)
    {
        x[i]=x[i-1]+deltax;
    }
    double sum=0; //initializing the sum of function values
    for(int i=0;i<N;i++)
    {
        if(i==0||i==N-1) sum+=f(x[i]);
        else sum+=2*f(x[i]);
    }
    double ans=(0.5)*deltax*sum; // the final result
    return ans;
}

double simpson(int N)
{
    double deltax=(double)(3-0)/(double)(N-1); //Dividing the region into N-1 equally spaced points
    vector<double> x(N);    //An array that includes the value of x at every node
    x[0]=0;
    for(int i=1;i<N;i++)
    {
        x[i]=x[i-1]+deltax;
    }
    double sum=0;
    for(int i=0;i<N;i++)
    {
        if(i==0||i==N-1) sum+=f(x[i]);
        else{
            if(i%2) sum+=4*f(x[i]);
            else sum+=2*f(x[i]);
        }
    }
    double ans=(1/3.00)*deltax*sum; // the final answer
    return ans;
}

int main()
{
    vector<int> N={5,11,41,101};

    double analytical= 1-4*exp(-3); // From Integration by Parts , exact solution = -(1+x)*e^(-x) + c ;

    printf("N  Exact  Trapezoidal Rule  Simpson's Rule\n");

    for(int i=0;i<4;i++)
    {
        printf("%d  ",N[i]);
        printf("%f  ",analytical);
        printf("%f  ",trapezoidal(N[i]));
        printf("%f\n",simpson(N[i]));
    }
}
