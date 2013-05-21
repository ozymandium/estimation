/*!
*   \file bryson.hpp
*
*   Implement Bryson's trick to approximate the process noise covariance, and 
*   optionally the discretized system matrix.
*
*  Copyright (c) 2013
*  \author Cofield, Robert
*
*/

#ifndef _BOOST_UBLAS_BRYSON_
#define _BOOST_UBLAS_BRYSON_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

/* backbone function */
UDmatrix calculate_C_bryson ( UDmatrix Ac,
                                    UDmatrix G,
                                    UDmatrix Qc,
                                    double Ts,
                                    uint n )
{
    UDmatrix GQG = ublas::prod( ublas::prod(G,Qc), ublas::trans(G) );
    UDmatrix Atrans = ublas::trans(Ac);

    UDmatrix S (2*n, 2*n);
    // S11 = Ac
    for (uint i=0; i<n; ++i)
        for (uint j=0; j<n; ++j)
            S[i,j] = Ac[i,j];
    // S12 = GQG
    for (uint i=n; i<2*n; ++i)
        for (uint j=n; j<2*n; ++i)
            S[i,j] = GQG[i,j-n];
    // S21 = 0
    for (uint i=n; i<2*n; ++i)
        for (uint j=0; j<n; ++i)
            S[i,j] = 0f;
    // S22 = Atrans
    for (uint i=n; i<2*n; ++i)
        for (uint j=n; j<2*n; ++i)
            S[i,j] = Atrans[i-n,j-n];

    UDmatrix C (2*n, 2*n);
    C = S*Ts;

    return C;
}

/* actual api functions */
namespace boost { namespace numeric { namespace ublas {

/* 
    Calculate only Qd
    return 1 if Ad from bryson isn't equivalent to expm(Ac*Ts)
*/
int bryson (UDmatrix Ac,
            UDmatrix G,
            UDmatrix Qc,
            double Ts,
            UDmatrix Qd )
{
    unsigned n = Ac.size1;
    // TODO check that Ac is square
    UDmatrix C = calculate_C_bryson(Ac,G,Qc,Ts,n);

    return 0;
}


}}}

#endif