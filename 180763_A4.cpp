#include<bits/stdc++.h>
using namespace std;
#define pi 3.141592653589
int n;
vector<vector<double>> a; double alpha=1,lambda=9,L=1;

double deltat,deltax,Fo;
vector<double> b(n,0),x(n,0);

bool conv(vector<double> y, vector<double> b)  /// Convergence Criteria for Gauss Siedel ( norm(Ax-B) < epsilon )
{
    double sum=0;
    for(int i=0;i<n;i++)
    {
        double sum1=0;
        for(int j=0;j<n;j++)
        {
            sum1+=(a[i][j]*y[j]);
        }
        sum1-=b[i];
        sum+=sum1*sum1;
    }
    sum=sqrt(sum);
    if(sum<0.000001) return true; else return false;
}

double anal(double x,double t)  /// Function for the analytical solution
{
    if(x==0) return 1;
    double sum=0;
    for(int k=1;k<=100000;k++)   /// We have taken no. of terms here 100000 for better accuracy. I have attached the image file with 100 terms result
    {
        sum+=((2*k*pi)*sin(k*pi*x)/(lambda+(pi*pi*k*k)))*(1-exp(-t*(lambda+(pi*pi*k*k))));
    }
    return sum;
}

double d_anal(double x,double t) /// Function for gradient for analytical solution
{
    double e=L/(double)(n-1);
    return (anal(x+e,t)-anal(x,t))/e;   ///Forward difference formula at the left end
}

solve()
{
    deltax=L/(double)(n-1);
    Fo=alpha*deltat/(deltax*deltax);  /// Grid Fourier Number

    a.resize(0);
    a.resize(n,vector<double>(n,0));        /// We now define the coefficient matrix
    a[0][0]=1; a[n-1][n-1]=1;

    for(int i=1;i<n-1;i++)
    {
        a[i][i-1]=Fo; a[i][i]=-(2*Fo+lambda*deltat+1); a[i][i+1]=Fo;
    }


    b.resize(0);x.resize(0);
    b.resize(n,0);x.resize(n,0);

    b[0]=1;


    for(int q=1;q<=20;q++){
        while(1)    /// Loop for Gauss Siedel
        {
            vector<double> xnew(n);
            for(int i=0;i<n;i++)
            {
                double sum=0;
                for(int j=0;j<n;j++)
                {
                    if(j<i) sum+=a[i][j]*xnew[j];
                    else if(j>i) sum+=a[i][j]*x[j];
                }
                xnew[i]=(b[i]-sum)/a[i][i];
            }
            double sum1=0;
            for(int i=0;i<n;i++) sum1+=pow((xnew[i]-x[i]),2); /// Calculating the norm of the sum
            sum1=sqrt(sum1);

            if(sum1<0.00001 && conv(xnew,b)) break;
            else{
                x=xnew;
            }
        }
        for(int i=1;i<n-1;i++) b[i]=-x[i];

        printf("%d        \t%f\t%f\t %f \t\t  ",q,q*deltat,((x[1]-x[0])/deltax),d_anal(0,q*deltat));

        printf("%f\t\t %f\n",x[(n-1)/2],anal(0.5,q*deltat));

    }

}

int main()
{
    cout<<"n=101\t time step chosen=0.000535 sec\n";
    cout<<"Iteration Time Elapsed Gradient(Numerical) Gradient(Analytical) Midpoint Value(Numerical) Midpoint Value(Analytical)\n";
    n=101; deltat=0.000535; solve(); cout<<"\n";
    cout<<"n=201\t time step chosen=0.00054 sec\n";
    cout<<"Iteration Time Elapsed Gradient(Numerical) Gradient(Analytical) Midpoint Value(Numerical) Midpoint Value(Analytical)\n";
    n=201; deltat=0.00054; solve();
}
