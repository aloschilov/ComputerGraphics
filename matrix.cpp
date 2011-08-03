///////////////////////////////////////////////////////////////////////////////
/// @file Matrix.cpp
/// @brief This file contains implementation of Matrix class
///
/// Related Files: The definitions are in Matrix.h
/// Other related files are: types.h
///
/// Maintained by: NSTU
///
/// Copyright c 2007 aloschil. All right reserved.
///
/// CONFIDENTIALLY AND LIMITED USE
///
/// This software, including any software of third parties embodied herein,
/// contains information and concepts wich are confidential to Alexander
/// Loschilov and such third parties. This software is licensed for use
/// solely in accordance with the terms and conditions of the applicable
/// license agreement with Alexander Loschilov or his authorized distributor.
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

#include "matrix.h"

namespace Aloschil
{
const double EPSILON = 1e-10;
double Matrix::epsilon = EPSILON;
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
bool operator == (const Matrix &m1,const Matrix &m2)
{
        if(m1.x!=m2.x || m1.y!=m2.y)
        {
                return false;
        }
        for(unsigned int i=0;i<m1.x;++i)
                for(unsigned int j=0;j<m1.y;++j)
                if(m1.at(j,i)!=m2.at(j,i))
                        return false;
        return true;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix augment(const Matrix &v,const Matrix &omega)
{
        if(v.rows()!=omega.rows()) throw ErrorSize();

        Matrix resultMatrix;
        resultMatrix.setSize(v.rows(),v.columns()+omega.columns());
        for(unsigned int i=0;i<v.rows();++i)
                for(unsigned int j=0;j<v.columns();++j)
                {
                        resultMatrix.at(i,j) = v.at(i,j);
                }

        for(unsigned int i=0;i<v.rows();++i)
                for(unsigned int j=0;j<omega.columns();++j)
                {
                        resultMatrix.at(i,j+v.columns()) = omega.at(i,j);
                }
        return resultMatrix;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Vector(double x,double y)
{
        Matrix tmp;
        tmp.setSize(2,1);
        tmp[0]=x;
        tmp[1]=y;
        return tmp;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Vector4(double x,double y,double z,double a)
{
	Matrix tmp;
	tmp.setSize(4,1);
	tmp[0]=x;
	tmp[1]=y;
	tmp[2]=z;
	tmp[3]=a;
	return tmp;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix::Matrix(double _epsilon)
{
setRangeOfEquivalence(_epsilon);
 changed = true;
 field = 0;
 x = 0;
 y = 0;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix::Matrix(u_int xNew,u_int yNew,double *array,double _epsilon)

{
setRangeOfEquivalence(_epsilon);
 if(!xNew || !yNew)
 {
  field = 0;
  x = 0;
  y = 0;
 }
 else
 {
  x=xNew;
  y=yNew;
  field = new double [x*y];
  if(!field) throw ErrorNoMemory();
   for(u_int i=0;i<x*y;i++)
   field[i]=array[i];
 }
 changed = true;
 name = "Matrix";
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix::Matrix(const Matrix &M)
{
 changed=M.changed;
 determinant=M.determinant;
 x=M.x;
 y=M.y;
 name=M.name;
 epsilon=M.epsilon;
 field = new double [x*y];
 if(!field) throw ErrorNoMemory();
 for(u_int i=0;i<x*y;i++)
  field[i]=M.field[i];
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
double & Matrix::operator [](unsigned int LinePos)
{
 if(!field)
  throw ErrorNoInitBase();
 if(LinePos>x*y-1) throw  ErrorUsingOperator();
 return field[LinePos];
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
double & Matrix::at(u_int row,u_int column) const
{
 if(!field) throw ErrorNoInitBase();
 if(row*x + column > x*y-1) throw  ErrorUsingOperator();
 return field[row*x + column];
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator = (const Matrix &M)
{
 if(this == &M) return *this;

 DeleteMatrix();
 changed=M.changed;
 determinant=M.determinant;
 x=M.x;
 y=M.y;
 epsilon=M.epsilon;

 field = new double [x*y];
 if(!field) throw ErrorNoMemory();
 for(u_int i=0;i<x*y;++i)
  field[i]=M.field[i];
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator += (const Matrix &M)
{
 if(!field)
  throw ErrorNoInitBase();
 if(x!=M.x || y!=M.y)
  throw ErrorSizeDiff();
 for(u_int i=0;i<x*y;i++) field[i]+=M.field[i];

 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator *= (const Matrix &M)
{
 if(!field)
  throw ErrorNoInitBase();
 if(!M.field)
  throw ErrorNoInitSecond();
 if(x!=M.y)
  throw ErrorSizeDiff();

 double *c = new double[M.x*y];
 if(!c) throw ErrorNoMemory();
 for(u_int l=0;l<M.x*y;++l) c[l]=0;

 for(u_int i=0;i<y;i++)
  for(u_int j=0;j<M.x;j++)
   for(u_int n=0;n<x;n++)
    c[i*M.x+j]+=field[i*x+n]*M.field[M.x*n+j];
 x=M.x;

 delete [] field;
 field = c;
 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator -= (const Matrix &M)
{
 *this+=-M;
 return *this;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator /= (const Matrix &M)
{
 Matrix P(M);
 (*this)*= !P;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix  Matrix::operator + (const Matrix &M) const
{
 Matrix P= *this;
 P+= M;
 return P;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix  Matrix::operator * (const Matrix &M) const
{
 Matrix P = *this;
 P*=M;
 return P;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix  Matrix::operator - (const Matrix &M) const
{
 Matrix P = *this;
 P-=M;
 return P;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix  Matrix::operator / (const Matrix &M) const
{
 Matrix P = *this;
 P/=M;
 return P;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix &Matrix::operator ! ()
{
 //     ���� ��� ����� ���.
    if(!field)
        throw ErrorNoInitBase();
    if(x!=y)
        throw ErrorSizeSqr();
    bool result;
    double *b;
    double *c;
    double d;

    u_int i;
    u_int j;
    u_int k;
    u_int p;

    double w;
    double y;
    double *z;

    d = 1;
    result = true;
    u_int n = x;

    b = new double [n];
    if(!b) throw ErrorNoMemory();
    for(i=0;i<n;i++)b[i]=0;

    c = new double [n];
    if(!c) throw ErrorNoMemory();
    for(i=0;i<n;i++)c[i]=0;

    z = new double [n];
    if(!z) throw ErrorNoMemory();
    for(i=0;i<n;i++)z[i]=0;
    for(i = 1; i <= n; i++)
    {
        z[i-1] = i;
    }
    i = 1;
    do
    {
        k = i;
        y = field[(i-1)*x + i - 1];
        if( i+1<=n )
        {
            j = i+1;
            do
            {
                w = field[(i-1)*x + j - 1];
                if( fabs(w)>fabs(y) )
                {
                    k = j;
                    y = w;
                }
                j = j+1;
            }
            while(j<=n);
        }
        d = d*y;
        if( fabs(y)<=epsilon )
        {
            result = false;
            i = n+1;
        }
        else
        {
            y = double(1)/double(y);
            for(j = 1; j <= n; j++)
            {
                c[j-1]=field[(j-1)*x + k - 1];
                field[(j-1)*x + k - 1] = field[(j-1)*x + i - 1];
                field[(j-1)*x + i - 1]= - c[j-1]*y;
                b[j-1]=field[(i-1)*x + j - 1]*y;
                field[(i-1)*x + j - 1] = b[j-1];
            }
            j = (u_int) z[i-1];
            z[i-1] = z[k-1];
            z[k-1] = j;
            field[(i-1)*x + i - 1] = y;
            for(k = 1; k <= n; k++)
            {
                if( k!=i )
                {
                    for(j = 1; j <= n; j++)
                    {
                        if( j!=i )
                        {
         field[(k-1)*x + j - 1]=field[(k-1)*x + j - 1] - b[j-1]*c[k-1];
                        }
                    }
                }
            }
        }
        i = i+1;
    }
    while(i<=n);
    if( result )
    {
        for(i = 1; i <= n; i++)
        {
            k = (u_int) z[i-1];
            while(k!=i)
            {
                for(j = 1; j <= n; j++)
                {
                    w = field[(i-1)*x + j - 1];
                    field[(i-1)*x + j - 1] = field[(k-1)*x + j - 1];
                    field[(k-1)*x + j - 1] = w;
                }
                p = (u_int) z[i-1];
                z[i-1] = z[k-1];
                z[k-1] = p;
                d = -d;
                k = (u_int) z[i-1];
            }
        }
    }
    if(b) delete b;
    if(c) delete c;
    if(z) delete z;
    if(!result) throw ErrorMath();
    changed = true;
    return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator = (const double &m)
{
 DeleteMatrix();
 x=1;
 y=1;
 field = new double [1];
 if(!field) throw ErrorNoMemory();
 field[0]=m;
 changed = true;
  return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator += (const double &m)
{
 if(!field)
  throw ErrorNoInitBase();
 for(u_int i=0;i<x*y;i++) field[i]+=m;
 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator *= (const double &m)
{
 if(!field)
  throw ErrorNoInitBase();
 for(u_int i=0;i<x*y;i++) field[i]*=m;
 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator -= (const double &m)
{
 if(!field)
  throw ErrorNoInitBase();
 for(u_int i=0;i<x*y;i++) field[i]*=m;
 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix & Matrix::operator /= (const double &m)
{
 if(!field)
  throw ErrorNoInitBase();
 for(u_int i=0;i<x*y;i++) field[i]/=m;
 changed = true;
 return (*this);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Matrix::operator + (const double &m)
{
 Matrix C = *this;
 C+=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Matrix::operator * (const double &m)
{
 Matrix C = *this;
 C*=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Matrix::operator - (const double &m)
{
 Matrix C = *this;
 C-=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Matrix::operator / (const double &m)
{
 Matrix C = *this;
 C/=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix operator + (const double &m,const Matrix &M)
{
 Matrix C = M;
 C+=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix operator * (const double &m,const Matrix &M)
{
 Matrix C = M;
 C*=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix operator - (const double &m,const Matrix &M)
{
 Matrix C = M;
 C+=-m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix operator / (const double &m,const Matrix &M)
{
 Matrix C = M;
 !C;
 C*=m;
 return C;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::DeleteMatrix()
{
 if(field) delete[] field;
 x = 0;
 y = 0;
 field = 0;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
double Matrix::getDeterminant()
{
 if(!field)
        throw ErrorNoInitBase();
 if(x != y)
        throw ErrorSizeSqr();
 u_int col,row,i;
 double  x1,x2,t,t1,t2,c,s,denomin;
 if(!changed) return determinant;
 for(col = 0; col < x; col++)
    {
        for(row = col+1; row < y; row++)
        {
            x1 = field[col*x + col];
            x2 = field[row*x + col];
            if( x1==0&&x2==0 )
            {
                continue;
            }
            if( fabs(x1)>fabs(x2) )
            {
                t = x2/x1;
                denomin = fabs(x1)*sqrt(1+t*t);
            }
            else
            {
                t = x1/x2;
                denomin = fabs(x2)*sqrt(1+t*t);
            }
            c = x1/denomin;
            s = -x2/denomin;
            for(i = 0; i < x; i++)
            {
                t1 = field[col*x + i];
                t2 = field[row*x + i];
                field[col*x + i] = c*t1-s*t2;
                field[row*x + i] = s*t1+c*t2;
            }
        }
    }
    determinant = 1;
    for(i = 0; i < x; i++)
    {
        determinant*=field[i*x + i];
    }
    changed=false;
    return determinant;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::print()
{
 if(!field) throw ErrorNoInitBase();
 for(u_int i=0;i<y;i++)
 {
  for(u_int j=0;j<x;j++)
  {
   printf("|%10.3g",field[i*x+j]);
  }
  putchar('\n');
 }
 putchar('\n');
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix::~Matrix()
{
 DeleteMatrix();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix Matrix::operator - () const
{
 if(!field) throw ErrorNoInitBase();
 Matrix M = (*this);
 M*=-1;
 return M;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::InsertMatrix(const Matrix &v,int i,int j)
{
 for(int k = i;k<(int)y && k-i<(int)v.y;++k)
  for(int l = j;l<(int)x && l-j<(int)v.x;++l)
  if(k>=0 && l>=0)
   field[k*x + l] = v.field[(k-i)*v.x + (l-j)];
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::setSize(u_int rows,u_int columns)
{
 if(!rows && !columns) {if(field) delete [] field;field=0;x=y=0;return;}

 if(rows!=y || columns!=x)
 {
  u_int y_m = (rows<y) ? rows : y;
  u_int x_m = (columns<x) ? columns : x;

  double *new_field = new double[rows*columns];

  if(!new_field) throw ErrorNoMemory();
  for(u_int i1=0;i1<rows;++i1)
   for(u_int j1=0;j1<columns;++j1)
   if(i1<y_m && j1<x_m && field)
    new_field[i1*columns + j1] = field[i1*x + j1];
   else
    new_field[i1*columns + j1] = 0;

  x = columns; y = rows;
  if(field) delete [] field;
  field = new_field;
 }
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
/*Matrix Matrix::Gauss()
{
 Matrix temp=*this;
 temp.Gauss_M();

 vector<double> ret;
 double c=0;
 for(u_int i=0;i<y;++i)
  ret.push_back(c=temp.field[x*i+(x-1)]);
 return ret;

} */
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::Gauss_M()
{
 if(!field)
        throw ErrorNoInitBase();
 if(x != y+1)
        throw ErrorSize();

 for(u_int j=0;j<y;++j)
   SwapLines(FindMaxAbsInColon(j,j),j);
 for(u_int k=0;k<y;++k)
 {
   if(fabs(at(k,k))< epsilon)
   {
    u_int i;
     for(i=k+1;i<y;++i)
     if(fabs(at(i,k))< epsilon)
     {
      SwapLines(i,k);
      break;
     }
     /*�� �������� ���*/
     if(i==y) throw ErrorMath();
   }

    for(u_int i=k+1;i<y;++i)
    {
     double c = at(i,k)/at(k,k);
     at(i,k) = 0;
     for(u_int j=k+1;j<x;++j)
      at(i,j)-=c*at(k,j);
    }
 }
 if(fabs(at(y-1,y-1))<epsilon) throw ErrorMath();
 for(int i=y-1;i>=0;--i)
 {
  for(u_int j=i+1;j<y;++j)
   at(i,x-1)-=at(i,j)*at(j,x-1);
  at(i,x-1)/=at(i,i);
 }
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
/*u_int Matrix::FindMaxInColon(u_int &index,u_int BeginingIndex)
{
 if(!field)
        throw ErrorNoInitBase();
 if(index>x-1 || BeginingIndex>y-1)
        throw ErrorUsingOperator();
 double Max=field[x*BeginingIndex + index];
 u_int IMax = BeginingIndex;
 for(u_int j=BeginingIndex;j<y;++j)
        if(field[x*j + index]>Max) {Max=field[x*j + index]; IMax = j;}
 return IMax;
} */
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
Matrix::u_int Matrix::FindMaxAbsInColon(u_int &index,u_int BeginingIndex)
{
 if(!field)
        throw ErrorNoInitBase();
 if(index>x-1 || BeginingIndex>y-1)
        throw ErrorUsingOperator();
 double Max=fabs(field[x*BeginingIndex + index]);
 u_int IMax = BeginingIndex;
 for(u_int j=BeginingIndex;j<y;++j)
        if(fabs(field[x*j + index])>Max) {Max=fabs(field[x*j + index]); IMax = j;}
 return IMax;
} 
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
/*u_int Matrix::FindMinInColon(u_int &index,u_int BeginingIndex)
{
	if(!field)
		throw ErrorNoInitBase();
	if(index>x-1 || BeginingIndex>y-1)
		throw ErrorUsingOperator();
	double Min=field[x*BeginingIndex + index];
	u_int IMin = BeginingIndex;
	for(u_int j=BeginingIndex;j<y;++j)
		if(field[x*j + index]<Min) {Min=field[x*j + index]; IMin = j;}
	return IMin;
} */
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
/*u_int Matrix::FindMinAbsInColon(u_int &index,u_int BeginingIndex)
{
	if(!field)
		throw ErrorNoInitBase();
	if(index>x-1 || BeginingIndex>y-1)
		throw ErrorUsingOperator();
	double Min=fabs(field[x*BeginingIndex + index]);
	u_int IMin = BeginingIndex;
	for(u_int j=BeginingIndex;j<y;++j)
		if(fabs(field[x*j + index])<Min) {Min=fabs(field[x*j + index]); IMin = j;}
	return IMin;
}*/
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::SwapLines(u_int _first,u_int _second)
{
 if(!field)
        throw ErrorNoInitBase();
 if(_first>y-1 || _second>y-1)
        throw ErrorUsingOperator();
 for(u_int j=0;j<x;++j)
 {
  double temp = field[x*_first + j];
  field[x*_first + j] = field[x*_second + j];
  field[x*_second + j] = temp;
 }
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
void Matrix::PrintToFile(const char *fileName)
{
        FILE *file;
        if((file=fopen(fileName,"w"))==NULL)
                return;
         if(!field) throw ErrorNoInitBase();
        for(u_int i=0;i<y;i++)
 {
  for(u_int j=0;j<x;j++)
  {
   fprintf(file,"|%10.3g",field[i*x+j]);
  }
  fprintf(file,"\n");
 }
 fprintf(file,"\n");
 fclose(file);

}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------

Matrix Matrix::getTransposed() const
{
        if(!field)
         throw ErrorNoInitBase();

        Matrix resultMatrix;
        resultMatrix.setSize(this->columns(),this->rows());

        for(unsigned int i=0;i<this->rows();++i)
                for(unsigned int j=0;j<this->columns();++j)
                {
                        resultMatrix.at(j,i) = this->at(i,j);
                }
        return resultMatrix;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to matrix.h for details.
//----------------------------------------------------------------------------
bool Matrix::isInitialized()
{
	if(!field) return false;
	return true;
}

} // namespace Aloschil



