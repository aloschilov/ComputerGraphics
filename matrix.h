#ifndef _Aloschil_Matrix_h_
#define _Aloschil_Matrix_h_

///////////////////////////////////////////////////////////////////////////////
/// @file Matrix.h
/// @brief This file contains declaration of Matrix class
///
/// Related Files: The implementation is in Matrix.cpp
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

#include <list>

namespace Aloschil
{
/// Detemines equivalence range of double value
extern const double EPSILON;

///////////////////////////////////////////////////////////////////////////////
/// @class MatrixError
/// @brief Base class for processing Matrix class exceptions
///
/// This class will be used while catching exceptions
/// situations within calculating block that is using
/// Matrix class in order to define alternative way of
/// solving a task.
///////////////////////////////////////////////////////////////////////////////
class MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorNoMemory
/// @brief Impossibility to allocate memory in method of Matrix
///////////////////////////////////////////////////////////////////////////////
class ErrorNoMemory: public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorNoInit
/// @brief Working with uninitilized instance of Matrix
///////////////////////////////////////////////////////////////////////////////
class ErrorNoInit : public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorNoInitBase
/// @brief Calling Matrix instance is uninitilized
///////////////////////////////////////////////////////////////////////////////
	class ErrorNoInitBase : public ErrorNoInit{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorNoInitSecond
/// @brief Parameter Matrix instance is uninitilized
///////////////////////////////////////////////////////////////////////////////
	class ErrorNoInitSecond : public ErrorNoInit{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorUsingOperator
/// @brief Uncorrect parameter is passed to method dealing with Matrix
///////////////////////////////////////////////////////////////////////////////
class ErrorUsingOperator : public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorSize
/// @brief Impossible to continue calc. because of the incompatible sizes
///////////////////////////////////////////////////////////////////////////////
class ErrorSize : public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorSizeSqr
/// @brief Square matrix is needed for operation, nut is not given
///////////////////////////////////////////////////////////////////////////////
	class ErrorSizeSqr : public ErrorSize{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorSizeDiff
/// @brief Operands are to be with equal sizes but are not
///////////////////////////////////////////////////////////////////////////////
	class ErrorSizeDiff : public ErrorSize{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorMath
/// @brief Operands are to be with equal sizes but are not
///////////////////////////////////////////////////////////////////////////////
class ErrorMath : public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorMathZeroVector
/// @brief Operands are to be non-zero vertors
///////////////////////////////////////////////////////////////////////////////
        class ErrorMathZeroVector : public ErrorMath{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorOpeningFile
/// @brief Try to initiate Matrix instance from file. Just unable to open
///////////////////////////////////////////////////////////////////////////////
class ErrorOpeningFile : public MatrixError{};
///////////////////////////////////////////////////////////////////////////////
/// @class ErrorSavingToFile
/// @brief Try to save Matrix instance to file and just unable to.
///////////////////////////////////////////////////////////////////////////////
class ErrorSavingToFile : public MatrixError{};

///////////////////////////////////////////////////////////////////////////////
/// @class Matrix
/// @brief This class deals with basic matrix operations
///////////////////////////////////////////////////////////////////////////////
class Matrix
{
////////////////////////////////////////////////////////////////////
/// @var u_int
///
/// @brief Just simplification of unsigned int.
////////////////////////////////////////////////////////////////////
                typedef unsigned int u_int;

                friend Matrix operator + (const double &,const Matrix &);
                friend Matrix operator * (const double &,const Matrix &);
                friend Matrix operator - (const double &,const Matrix &);
                friend Matrix operator / (const double &,const Matrix &);
                friend bool operator == (const Matrix&,const Matrix &);
                friend void Swap_d(Matrix&,Matrix&);
                friend void printList(const char *fileName,std::list<Matrix> listToPrint);
                friend void printListOfLists(const char *fileName,std::list<std::list<Matrix> >listToPrint);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Function to initialize matrix with column vector
///
/// Example:
/// <PRE>
/// Matrix A = Vector(5,7);
/// A.print();
/// </PRE>
/// Output result:<br>
/// <B><I> |5<br>|7 </I></B>
/// @param x � input � value to initialize first element of the vector.
/// @param y � input � value to initialize second element of the vector.
///
/// @returns � 2D Column vector represented as Matrix.
///
//////////////////////////////////////////////////////////////////////////////
                friend Matrix Vector(double x,double y);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Function to initialize matrix with column vector
///
/// Example:
/// <PRE>
/// Matrix A = Vector(5,7);
/// A.print();
/// </PRE>
/// Output result:<br>
/// <B><I> |5<br>|7 </I></B>
/// @param x � input � value to initialize first element of the vector.
/// @param y � input � value to initialize second element of the vector.
///
/// @returns � 2D Column vector represented as Matrix.
///
//////////////////////////////////////////////////////////////////////////////
                friend Matrix Vector4(double x,double y,double z,double a=1.0);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Function joins arrays side by side
///
/// Example:
/// <PRE>
/// Matrix v;
/// v.setSize(2,2);
/// // 1 2
/// // 3 4
/// v.at(0,0)=1;v.at(0,1)=2;
/// v.at(1,0)=3;v.at(1,1)=4;
///
/// Matrix omega;
/// omega.setSize(2,3);
/// // 0  10 20
/// // 30 40 50
/// omega.at(0,0)=0;  omega.at(0,1)=10; omega.at(0,2)=20;
/// omega.at(1,0)=30; omega.at(1,1)=40; omega.at(1,2)=50;
///
/// Matrix result = augment(v,omega);
/// // 1 2  0 10 20
/// // 3 4 30 40 50
///
/// </PRE>
/// @param v � input � First matrix to combine.
/// @param omega � input � Second matrix to combine.
/// @warning � Matrixes are to be with the same number of rows otherwise
/// ErrorSize exception is thrown.
///
/// @returns � Joined matrix.
///
//////////////////////////////////////////////////////////////////////////////
                friend Matrix augment(const Matrix &v,const Matrix &omega);
         public:
//////////////////////////////////////////////////////////////////////////////
/// @brief Constructor that creates uninitialized matrix.
///
/// No operations are possible with matrix created this way.
///
/// @param _epsilon � input � This value used to set range of equivalence.
//////////////////////////////////////////////////////////////////////////////
                Matrix( double _epsilon=EPSILON);
//////////////////////////////////////////////////////////////////////////////
/// @brief Constructor that creates matrix initialized with array.
///
/// No operations are possible with matrix created this way.
///
/// @param numberOfColumns � input � Number of columns in new matrix.
/// @param numberOfRows � input � Number of rows in new matrix.
/// @param pArray � input � Pointer to first element to linear array.
/// @param _epsilon � input � This value used to set range of equivalence.
/// @warning � All the numberOfColumns*numberOfRows elements are to be
/// initialized in array provided. In case of non-allocating memory
/// segmentation fault error is possible.
//////////////////////////////////////////////////////////////////////////////
                Matrix( u_int numberOfColumns,
                        u_int numberOfRows,
                        double *pArray,
                        double _epsilon=EPSILON);
//////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
///
/// @param matrixToCopy � input � Matrix to duplicate with.
//////////////////////////////////////////////////////////////////////////////
                Matrix(const Matrix& matrixToCopy);
//////////////////////////////////////////////////////////////////////////////
/// @brief Destructor
//////////////////////////////////////////////////////////////////////////////
                ~Matrix();

//////////////////////////////////////////////////////////////////////////////
/// @brief Destructor
//////////////////////////////////////////////////////////////////////////////

                double & operator [](unsigned int) ;

                Matrix & operator  = (const Matrix &);

                Matrix & operator += (const Matrix &);
                Matrix & operator *= (const Matrix &);
                Matrix & operator -= (const Matrix &);
                Matrix & operator /= (const Matrix &);

                Matrix  operator + (const Matrix &) const;
                Matrix  operator * (const Matrix &) const;
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Subtraction operator
///
/// Subtracts Matrix standing at the right side of the minus
/// from Matrix standing at the left side of the minus by substracting
/// corresponding elements
/// <br>Example:
/// <PRE>
/// Matrix A;
/// A.setSize(2,2);
/// // 1 2
/// // 3 4
/// A.at(0,0)=1;A.at(0,1)=2;
/// A.at(1,0)=3;A.at(1,1)=4;
/// Matrix B;
/// B.setSize(2,2);
/// // 0  10
/// // 30 40
/// B.at(0,0)=0; B.at(0,1)=10;
/// B.at(1,0)=30; B.at(1,1)=40;
///
/// Matrix result = A - B;
/// result.print();
/// </PRE>
/// Output result:<br>
/// <B><I> |1  |-8<br>|-27 |-36 </I></B>
///
/// @param rightSideMatrix � input � Matrix standing at the right side of sign.
/// @warning � Matrixes are to be with the same number of rows and columns
/// otherwise ErrorSizeDiff exception is thrown.
///
/// @returns � Result Matrix of the operation.
///
//////////////////////////////////////////////////////////////////////////////
                Matrix  operator - (const Matrix &rightSideMatrix) const;
                Matrix  operator / (const Matrix &) const;

                Matrix & operator = (const double &);

                Matrix & operator += (const double &);
                Matrix & operator *= (const double &);
                Matrix & operator -= (const double &);
                Matrix & operator /= (const double &);

                Matrix operator + (const double &);
                Matrix operator * (const double &);
                Matrix operator - (const double &);
                Matrix operator / (const double &);

                // inversion
                Matrix & operator !();
                // determinant
                double getDeterminant ();
                // mines operator
                Matrix operator - () const;

                double & at(u_int row,u_int column) const;

                Matrix getTransposed() const;

                                Matrix Gauss();

                u_int FindMaxInColon(u_int &index,u_int BeginingIndex = 0);
	        u_int FindMaxAbsInColon(u_int &index,u_int BeginingIndex = 0);

	        u_int FindMinInColon(u_int &index,u_int BeginingIndex = 0);
	        u_int FindMinAbsInColon(u_int &index,u_int BeginingIndex = 0);

                void SwapLines(u_int ,u_int );

                void InsertMatrix(const Matrix &,int,int);
                void setSize(u_int rows,u_int columns);


                const u_int& columns() const {return x;}
                const u_int& rows() const {return y;}

                void LoadFromFile(const char *);
                void PrintToFile(const char *);
                void print();

                char* GetName(){return name;}
                void SetName( char *str) {name = str;}

                static void setRangeOfEquivalence(const double &_epsilon)
                {
                        epsilon =_epsilon;
                }


		bool isInitialized();
        private:
                void DeleteMatrix();
                void Gauss_M();
                // Pointer to first element of linearized matrix
                double *field;

                // Range of equivalence
                static double epsilon;

                // Matrix has changed since last detetminant computation
                bool changed;

                // Matrix's determinant
                double determinant;

                u_int x,y;

                // Name of matrix
                char *name;
};

Matrix Vector(double x,double y);
}
#endif // _Aloschil_Matrix_h_
