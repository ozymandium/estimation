/*
    Discrete-time Linear Quadratic Regulator Solver

    http://web.ics.purdue.edu/~zheng33/sources/DARE_code.txt

    Uses LAPACK 
*/
#include <mathutils/dare.hpp>

matrix<double> dlqr(matrix<double>A,matrix<double>B,matrix<double>R,matrix<double>Q)
{
    //A,Q,X:n*n, B:n*m, R:m*m

    int n=A.size1(),m=B.size2();
    matrix<double,column_major>P(n,n);
    matrix<double,column_major>temp1(m,n);
    matrix<double,column_major>temp2(m,m);

    P=dare(A,B,R,Q);                //Solution of DARE
    temp1=prod(trans(B),P);             //B'*P,m*n
    temp2=R+prod(temp1,B);          //R+(B'*P)*B,m*m
    temp1=prod(trans(temp2),temp1); //m*n
    
    return prod(temp1,A);
}
