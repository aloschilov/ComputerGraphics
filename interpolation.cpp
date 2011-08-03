//---------------------------------------------------------------------------


#pragma hdrstop

#include <math.h>

#include "interpolation.h"

using namespace std;

namespace Aloschil
{
//----------------------------------------------------------------------------
// This is a protected API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Spline::T(unsigned int n,double t)
{
        Matrix result;
        result.setSize(n,1);
        for(unsigned int i=0;i<n+1;++i)
        {
                result[i]=pow(double(i),t);
        }
        return result;
}
//----------------------------------------------------------------------------
// This is a protected API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Spline::V(unsigned int n,double k)
{
        Matrix result;
        result.setSize(n,1);
        result[0]=0;
        for(unsigned int i=1;i<n+1;++i)
        {
                result[i] =double(i) * pow(k,double(i-1));
        }
		return result;
}
//----------------------------------------------------------------------------
// This is a protected API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Spline::Q(unsigned int N)
{
        Matrix result;
        result.setSize(4*N,4*N);
        for(unsigned int i=0;i<N;++i)
        {
                for(unsigned int j=0;j<3;++j)
                {
                        result.at(i*4,i*3 + j) = 1;
                        for(unsigned int k=1;k<4;++k)
                        {
                                result.at(i*4 + k,i*3 + j) = pow(double(j),double(k));
                        }
                }
        }
        return result;
}
//----------------------------------------------------------------------------
// This is a protected API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Spline::initializeWith(vector<Matrix> U)
{
        Q(3).PrintToFile("Q3.mtx");
}


}
//---------------------------------------------------------------------------

#pragma package(smart_init)
