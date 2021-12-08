#include<bits/stdc++.h>
using namespace std;

vector<double> q(6),b(6);   //q is the vector of variables and b is the residual function vector

double e=1e-6; //Taking the value of epsilon for central difference formula

vector<vector<double>> a(6,vector<double>(6,0)); //The matrix for the Jacobian

double F1(){ return q[0]*q[0]+q[1]*q[1]+q[2]*q[2]-14; }       //Defining the 6 functions

double F2(){ return 2*q[1]*q[1]+q[2]*q[2]+2*q[3]*q[3]-35; }

double F3(){ return q[2]*q[2]+q[3]*q[3]-2*q[4]*q[4]-10; }

double F4(){ return q[0]+q[4]-3; }

double F5(){ return q[3]+q[5]-4; }

double F6(){ return q[1]+q[5]-3; }

void Constructa()   //Function to construct the Jacobian Matrix
{
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            vector<double> x=q; //Using a temporary vector to store the q vector
            if(i==0){
                q[j]=x[j]+e; //Using the Central Difference Formula
                a[i][j]=F1(); //F(x+e) --- Value of the Function at x+e
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F1(); // F(x-e) --- Value of Function at x-e
                a[i][j]/=2*e;   //Dividing by 2e
            }
            if(i==1){                   //Likewise, for other rows of the Jacobian
                q[j]=x[j]+e;
                a[i][j]=F2();
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F2();
                a[i][j]/=2*e;
            }
            if(i==2){
                q[j]=x[j]+e;
                a[i][j]=F3();
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F3();
                a[i][j]/=2*e;
            }
            if(i==3){
                q[j]=x[j]+e;
                a[i][j]=F4();
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F4();
                a[i][j]/=2*e;
            }
            if(i==4){
                q[j]=x[j]+e;
                a[i][j]=F5();
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F5();
                a[i][j]/=2*e;
            }
            if(i==5){
                q[j]=x[j]+e;
                a[i][j]=F6();
                q[j]=x[j]-e;
                a[i][j]=a[i][j]-F6();
                a[i][j]/=2*e;
            }
            q=x;
        }
    }
}

void Constructb()       // Constructing the residual Vector
{
    b[0]=-F1(); b[1]=-F2(); b[2]=-F3(); b[3]=-F4(); b[4]=-F5();
    b[5]=-F6();
}


void GaussElimination()
{
    int n=6;
    for(int k=0;k<n-1;k++)
    {
        for(int i=k+1;i<n;i++)              // Loop to convert the Jacobian to an upper triangular matrix // Forward Elimination
        {
            double factor=a[i][k]/a[k][k];
            for(int j=0;j<n;j++)
            {
                a[i][j]=a[i][j]-factor*a[k][j];
            }
            b[i]=b[i]-factor*b[k];
        }
    }
    vector<double> x(n);
    for(int i=n-1;i>=0;i--)             // Back Substitution
    {
        double sum=b[i];
        for(int j=i+1;j<n;j++)
        {
            sum-=a[i][j]*x[j];
        }
        x[i]=sum/a[i][i];
    }
    for(int i=0;i<n;i++) b[i]=x[i];
}

int main()
{
    for(int i=0;i<6;i++) q[i]=5.00; //Initializing the guess value for q


    printf("Iteration Index \tGuessed Value(Q1--Q6)       \t\t\t\tCorrections(deltaQ1--deltaQ6)     \t\t\t\tUpdates(Qi+deltaQi)\n");

    for(int i=0;i<9;i++)
    {
        Constructa();
        Constructb();
        vector<double> deltax(6); // The deltaQ vector
        GaussElimination();
        deltax=b;
        cout<<i+1<<"\t\t\t";
        for(int j=0;j<6;j++)
        {
            printf("%f,",q[j]);
            q[j]=q[j]+deltax[j];
        }
        cout<<"\t";
        for(int j=0;j<6;j++)
        {
            printf("%f,",deltax[j]);
        }
        cout<<"\t";
        for(int j=0;j<6;j++)
        {
            printf("%f,",q[j]);
        }
        cout<<"\n";
    }

}
