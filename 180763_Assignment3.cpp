#include<bits/stdc++.h>
using namespace std;

vector<vector<double>> a(5,vector<double>(5,0)); // Coefficient matrix
double F; // Grid Fourier Number

void GaussElimination()
{
    double temp;
    int order=5;
    vector<vector<double>>matrix(order, vector<double>(2*order,0));
    for(int i = 0; i < order; i++){
        for(int j = 0; j <order; j++){
            matrix[i][j] = a[i][j];
        }
    }
    // Create the augmented matrix
    // Add the identity matrix
    // of order at the end of original matrix.
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < 2 * order; j++) {

            // Add '1' at the diagonal places of
            // the matrix to create a identity matirx
            if (j == (i + order))
                matrix[i][j] = 1;
        }
    }

    // Interchange the row of matrix,
    // interchanging of row will start from the last row
    for (int i = order - 1; i > 0; i--) {
        // Directly swapping the rows using pointers saves
        // time

        if (matrix[i - 1][0] < matrix[i][0]) {
            vector<double> tmp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = tmp;
        }
    }

    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 0; i < order; i++) {

        for (int j = 0; j < order; j++) {

            if (j != i) {

                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {

                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }

    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {

        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {

            matrix[i][j] = matrix[i][j] / temp;
        }
    }
    vector<vector<double>>res(order, vector<double>(order,0));
    for(int i = 0; i<order; i++){
        for(int j = 0; j < order; j++){
            res[i][j] = matrix[i][j+order];
        }
    }

    a=res;
}

void constructDDa() // Function to find the coefficient matrix for Dirichlet Dirichlet Condition
{
    a={{1,0,0,0,0},
    {F,-(1+2*F),F,0,0},
    {0,F,-(1+2*F),F,0},
    {0,0,F,-(1+2*F),F},
    {0,0,0,0,1}};

    for(int i=0;i<5;i++)     //Normalizing the matrix
    {
        double x=-1e9,y=-1e9;
        for(int j=0;j<5;j++)
        {
            if(x<abs(a[i][j]))   //Finding the element in a row with maximum absolute value
            {
                x=abs(a[i][j]); y=a[i][j];
            }
        }
        for(int j=0;j<5;j++)
        {
            a[i][j]=a[i][j]/y;  //Dividing the elements with maximum value in the row
        }
    }
}

void constructDNa() // Function to find the coefficient matrix for Dirichlet Neumann Condition
{
    a={{1,0,0,0,0},
    {F,-(1+2*F),F,0,0},
    {0,F,-(1+2*F),F,0},
    {0,0,F,-(1+2*F),F},
    {0,0,0,-1,1}};

    for(int i=0;i<5;i++)     //Normalizing the matrix
    {
        double x=-1e9,y=-1e9;
        for(int j=0;j<5;j++)
        {
            if(x<abs(a[i][j]))   //Finding the element in a row with maximum absolute value
            {
                x=abs(a[i][j]); y=a[i][j];
            }
        }
        for(int j=0;j<5;j++)
        {
            a[i][j]=a[i][j]/y;  //Dividing the elements with maximum value in the row
        }
    }
}

double Norm()
{
    double ret=-1e9;
    for(int i=0;i<5;i++)
    {
        double x=0;
        for(int j=0;j<5;j++)
        {
            x=x+abs(a[i][j]);
        }
        ret=max(ret,x);
    }
    return ret;
}



int main()
{
    vector<double> Fo={0.1,0.25};
    cout<<"Boundary Condition   \t  Fo  \t Norm(A) \t  Norm(A_inverse) \t  Condition Number\n";
    for(int i=0;i<2;i++)
    {
        cout<<"D-D  \t\t\t"; cout<<Fo[i]<<"\t";
        F=Fo[i];
        constructDDa();
        double n1,n2;
        printf("%.6f \t\t",Norm());
        n1=Norm();
        GaussElimination();
        cout<<Norm()<<" \t\t";
        n2=Norm();
        double CN=n1*n2;
        cout<<CN<<"\n";

        cout<<"D-N  \t\t\t"; cout<<Fo[i]<<"\t";
        constructDNa();
        printf("%.6f \t\t",Norm());
        n1=Norm();
        GaussElimination();
        cout<<Norm()<<" \t\t";
        n2=Norm();
        CN=n1*n2;
        cout<<CN<<"\n";
    }
}

