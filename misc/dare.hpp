/*
    Discrete-time Algebraic Riccati Equation solver

    http://web.ics.purdue.edu/~zheng33/sources/DARE_code.txt

    Uses LAPACK

*/           

matrix<double> dare(matrix<double> A, matrix<double> B, matrix<double> R,matrix<double> Q)
{
    //A,Q,X:n*n, B:n*m, R:m*m

    int n=A.size1(),m=B.size2();
    matrix<double,column_major>Z11(n,n);
    matrix<double,column_major>Z12(n,n);
    matrix<double,column_major>Z21(n,n);
    matrix<double,column_major>Z22(n,n);
    matrix<double,column_major>temp1(m,m);
    matrix<double,column_major>temp2(n,n);
    matrix<double,column_major>temp3(m,n);

    Z11=inv(A);                       //inv(A)
    temp1=inv(R);                     //inv(R)
    Z21=prod(Q,Z11);                  //Q*inv(A)
    temp2=prod(B,temp1,trans(B));     //B*inv(R)*B'
    Z12=prod(Z11,temp2);
    Z22=trans(A)+prod(Z21,temp2);

    //construct the Z with Z11,Z12,Z21,Z22
    matrix<double,column_major>Z(2*n,2*n);
    subrange(Z,0,n,0,n)=Z11;
    subrange(Z,n,2*n,0,n)=Z21;
    subrange(Z,0,n,n,2*n)=Z12;
    subrange(Z,n,2*n,n,2*n)=Z22;

    vector<std::complex<double> >eig(2*n);

    int n_Z=traits::matrix_size1(Z);
    matrix<double,column_major> VL(n_Z,n_Z);
    matrix<double,column_major> VR(n_Z,n_Z);

    lapack::geev(Z,eig,&VL,&VR,lapack::optimal_workspace());

    //Order the eigenvectors
    //move e-vectors correspnding to e-value outside the unite circle to the left
    matrix<double,column_major>U11(n,n);
    matrix<double,column_major>U21(n,n);
    matrix<double,column_major>tempZ(n_Z,n);

    int c1=0;

    for (int i=0;i<n_Z;i++) {
        if ((eig(i).real()*eig(i).real()+eig(i).imag()*eig(i).imag())>1) { //outside the unite cycle
            column(tempZ,c1)=column(VR,i);
            c1++;
        }
    }

    U11=subrange(tempZ,0,n,0,n);
    U21=subrange(tempZ,n,n_Z,0,n);

    return prod(U21,inv(U11));
}


