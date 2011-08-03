#ifndef _Aloschil_Linalg_h
#define _Aloschil_Linalg_h

///////////////////////////////////////////////////////////////////////////////
/// @file Matrix.h
/// @brief This file contains definition of Matrix class
///
/// Related Files: The implementation is in LinAlg.cpp
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

#include <vector>
#include <list>

#include "matrix.h"

namespace Aloschil
{
extern const double deg;

///////////////////////////////////////////////////////////////////////////////
/// @class NoCrossPoint
/// @brief Unable to get cross point of some elements
///////////////////////////////////////////////////////////////////////////////
class NoCrossPoint : public MatrixError{};

////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
/// @enum SectionsPositionalRelationship
/// @brief Characterizes positional relationship of two sections
///
/// These characteristics are directional.
///
/// <A NAME="SectionsPositionalRelationship">
///  <IMG src=images/SectionsPositionalRelationship_enum.jpg ALIGN=left>
/// </A>
///
/// On image (1) CD is TOUCHING, AB is TOUCHED.<br>
/// On image (2) CD is TOUCHED, AB is TOUCHING.<br>
/// On image (3) AB and CD are both CROSSING_OR_CROSSED.<br>
/// On image (4) AB and CD are both COLLINEAR.<br>
/// On image (5) AB and CD are both ADJACENT.<br>
/// On image (6) AB and CD are both NONINTERSECTING.
///
////////////////////////////////////////////////////////////////////
enum SectionsPositionalRelationship
{
        COLLINEAR,
        CROSSING_OR_CROSSED,
        TOUCHED,
        TOUCHING,
        ADJACENT,
        NONINTERSECTING
};

//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Computes angle in between 2D column vectors represented as Matrixes
///
/// The value returned is always between -M_PI and M_PI, including M_PI.
/// Returns 0 if one of vectors is zero-vector.
///
/// Example 1:
/// <PRE>
///        Matrix v;
///        v.setSize(2,1);
///        v[0]=1;
///        v[1]=3;
///
///        Matrix omega;
///        omega.setSize(2,1);
///        omega[0]=-1;
///        omega[1]=-0.000001;
///
///        printf("ang(v,omega) = %f\n",ang(v,omega)/deg);
///        // ang(v,omega) = 108.435006
///  </PRE>
///
/// @param v � input � Vector to count angle from.
/// @param omega � input � Vector to count angle to.
///
/// @returns � Angle in rand between vectors
///
//////////////////////////////////////////////////////////////////////////////
        double ang(const Matrix &v,const Matrix &omega);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates polar angle for 2D vector in Cartesian coordinate system
///
/// Calculates the angle (in radians) from the x-axis to a line containing the
/// origin (0,0) and the point (x,y). The value returned is always between
/// -M_PI and M_PI, including M_PI.
///
/// Example:
/// <PRE>
///        Matrix p;
///        p.setSize(2,1);
///        p[0]=-1;
///        p[1]=-0.000001;
///
///        printf("atan(p) = %f\n",atan(p)*deg);
///        // atan(p) = -179.999943
/// </PRE>
///
/// @param p � input � Vector to calculate polar angle of.
///
/// @returns � Returns the angle (in radians).
///
//////////////////////////////////////////////////////////////////////////////
        double atan(const Matrix &p);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates parameters (t - <A HREF="#tParameterImage">see image</A>) of cross point.
///
///  <A NAME="tParameterImage">
///  <IMG src=images/getCrossingParameters.jpg ALIGN=left>
///  </A>
///
/// Example 1:
/// <IMG src=images/getCrossingParameters_example_1.jpg ALIGN=left>
/// <PRE>
///        Matrix A = Vector(1,2);
///        Matrix B = Vector(5,4);
///        Matrix C = Vector(5,1);
///        Matrix D = Vector(2,4);
///
///        Matrix T = getCrossingParameters(A,B,C,D);
///        T.print();
/// </PRE>
/// Output 1:<B><I><br>| 0.5<br> | 0.667</I></B>
///
/// <br>Example 2:
/// <IMG src=images/getCrossingParameters_example_2.jpg ALIGN=left>
/// <PRE>
///        Matrix A = Vector(1,2);
///        Matrix B = Vector(5,4);
///        Matrix C = Vector(5,1);
///        Matrix D = Vector(4,2);
///
///        Matrix T = getCrossingParameters(A,B,C,D);
///        T.print();
/// </PRE>
/// Output 2:<B><I><br>| 0.5<br> | 2</I></B>
///
/// @param A � input � First point of first line segment.
/// @param B � input � Second point of first line segment.
/// @param C � input � First point of second line segment.
/// @param D � input � Second point of second line segment.
/// @warning � (B-A) and (C-D) vectors are to be canted otherwise
/// ErrorMath exception is thrown.
///
/// @returns � 2D column vector with t-parameters as coordinates.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix getCrossingParameters(const Matrix &A,const Matrix &B,
                                     const Matrix &C,const Matrix &D);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates cross point (<A HREF="#getCrosspoint">red one on image</A>) of two lines passing throught two points
///
///  <A NAME="getCrosspoint">
///  <IMG src=images/getCrosspoint.jpg ALIGN=left>
///  </A>
///
/// @param A � input � The first point the first line passing through.
/// @param B � input � The second point the first line passing through.
/// @param C � input � The first point the second line passing through.
/// @param D � input � The second point the second line passing through.
/// @warning � (B-A) and (C-D) vectors are to be canted otherwise
/// ErrorMath exception is thrown.
///
/// @returns � Crosspoint as 2D column vector.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix getCrossPoint(   const Matrix &A,const Matrix &B,
                                const Matrix &C,const Matrix &D);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates cross point (<A HREF="#getSectionsCrosspoint">red one on image</A>) of two lines passing throught two points
///
///  <A NAME="getSectionsCrosspoint">
///  <IMG src=images/getSectionsCrosspoint.jpg ALIGN=left>
///  </A>
///
/// @param A � input � The first point the first section.
/// @param B � input � The second point the first section.
/// @param C � input � The first point the second section.
/// @param D � input � The second point the second section.
/// @warning � (B-A) and (C-D) vectors are to be canted otherwise
/// ErrorMath exception is thrown.
/// If there is no cross point NoCrossPoint exception is thrown.
///
/// @returns � Crosspoint as 2D column vector.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix getSectionsCrossPoint(   const Matrix &A,const Matrix &B,
                                        const Matrix &C,const Matrix &D);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates cross point (<A HREF="#getCrosspoint">red one on image</A>) of two lines passing throught two points
///
///  <A NAME="getSectionsCrosspoint">
///  <IMG src=images/getSectionsCrosspoint.jpg ALIGN=left>
///  </A>
///
/// @param A � input � The first point the first section.
/// @param B � input � The second point the first section.
/// @param C � input � The first point the second section.
/// @param D � input � The second point the second section.
/// @warning � (B-A) and (C-D) vectors are to be canted otherwise
/// ErrorMath exception is thrown.
/// If there is no cross point NoCrossPoint exception is thrown.
///
/// @returns � Crosspoint as 2D column vector.
///
//////////////////////////////////////////////////////////////////////////////
        double getDistance(const Matrix &A, const Matrix &B);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Classifies section CD with respect to AB
///
/// @see SectionsPositionalRelationship for more information.
///
/// @param A � input � The first point the first section.
/// @param B � input � The second point the first section.
/// @param C � input � The first point the second section.
/// @param D � input � The second point the second section.
/// @param rangeOfEquivalence � input � Distance between points where they
/// are considered to be equal.
/// @param rangeOfCollinearity � input � Range of angle within which lines
/// are considered to be collinear.
///
/// @returns � Value that stands for class.
///
//////////////////////////////////////////////////////////////////////////////
        SectionsPositionalRelationship getPositionalRelationship(
                                const Matrix &A,const Matrix &B,
                                const Matrix &C,const Matrix &D,
                                double rangeOfEquivalence,
                                double rangeOfCollinearity);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Classifies section CD with respect to AB
///
/// @see SectionsPositionalRelationship for more information.
///
/// @param A � input � The first point the first section.
/// @param B � input � The second point the first section.
/// @param C � input � The first point the second section.
/// @param D � input � The second point the second section.
/// @param rangeOfEquivalence � input � Distance between points where they
/// are considered to be equal.
/// @param rangeOfCollinearity � input � Range of angle within which lines
/// are considered to be collinear.
///
/// @returns � Value that stands for class.
///
//////////////////////////////////////////////////////////////////////////////
        std::vector<Matrix> getRandomPolygon(   const Matrix &T,
                                                double a, double b,
                                                double teta);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Classifies section CD with respect to AB
///
/// @see SectionsPositionalRelationship for more information.
///
/// @param A � input � The first point the first section.
/// @param B � input � The second point the first section.
/// @param C � input � The first point the second section.
/// @param D � input � The second point the second section.
/// @param rangeOfEquivalence � input � Distance between points where they
/// are considered to be equal.
/// @param rangeOfCollinearity � input � Range of angle within which lines
/// are considered to be collinear.
///
/// @returns � Value that stands for class.
///
//////////////////////////////////////////////////////////////////////////////
        std::vector<Matrix> getRound(   const Matrix &a,
                                        const Matrix &b,
                                        const Matrix &c,
                                        double r,
                                        double rho,
                                        int numberOfSegments,
                                        double rangeOfCollinearity=EPSILON);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function generates set of points representing an arc
///
///
///
/// @param o � input � The first point the first section.
/// @param r � input � The second point the first section.
/// @param alfa � input � The first point the second section.
/// @param redoubledPhi � input � The second point the second section.
/// @param numberOfSegments � input � Distance between points where they
/// are considered to be equal.
///
/// @returns � Value that stands for class.
///
//////////////////////////////////////////////////////////////////////////////
        std::vector<Matrix> getArc(     const Matrix &o,
                                        double r,
                                        double alfa,
                                        double redoubledPhi,
                                        int numberOfSegments);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Used to rotate vector on 90 degrees counterclockwise.
///
/// @param vectorToRotate � input � 2D column vector.
///
/// @returns � 2D column vector rotated 90 degrees counterclockwise.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix Lrot(const Matrix &vectorToRotate);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Used to rotate vector on 90 degrees clockwise.
///
/// @param vectorToRotate � input � 2D column vector.
///
/// @returns � 2D column vector rotated 90 degrees clockwise.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix Rrot(const Matrix &vectorToRotate);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Used to get normalized vector.
///
/// @param vectorToNormalize � input � multidimensional column vector.
///
/// @returns � Normalized column vector with the same dimensionality as
/// for the input one.
///
//////////////////////////////////////////////////////////////////////////////
        Matrix normalize(const Matrix &vectorToNormalize);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Basic functions
///
/// @brief Used to get the sign of the parameter passed.
///
/// Example:
/// <PRE>
///        double a = -3423;
///        double b = 3423;
///        double c = 0;
///
///        printf("Value returned by sing(a) is %f.\r\n",sign(a));
///        // Value returned by sing(a) is -1.000000.
///        printf("Value returned by sing(b) is %f.\r\n",sign(b));
///        // Value returned by sing(b) is 1.000000.
///        printf("Value returned by sing(c) is %f.\r\n",sign(c));
///        // Value returned by sing(c) is 0.000000.
/// </PRE>
///
/// @param x � input � The variable to get the sign of.
///
/// @returns � Returns 0 if x = 0, 1 if x > 0, and 1 otherwise.
///
//////////////////////////////////////////////////////////////////////////////
        double sign(double x);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Compiles list of closed circuits from list of segments
///
/// Very strict constraints are imposed on segments circuits composed from.
/// It's to be no multiple segments.
/// It's to be no segments with deadlock vertexes.
/// It's to be only one paar on equivalent points.
///
/// @param listOfSegments - input - list of segments to extract circuits from
/// @param rangeOfEquivalence � input � Distance between points where they
/// are considered to be equal.
///
/// @returns - Returns list of closed circuits if all the conditions specified
/// above are met and empty list otherwise.
///
//////////////////////////////////////////////////////////////////////////////
        std::vector<std::vector<Matrix> > getClosedCircuits(
                                std::vector<std::vector<Matrix> > listOfSegments,
                                double rangeOfEquivalence=EPSILON);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Calculates cross product of two three-element column vectors
///
/// Direction of the cross-product is orthogonal to u and v in the direction
/// determined by the right hand rule
///
/// @param u - input - first three-element column vector to calculate cross
/// product of
/// @param v � input � second three-element column vector to calculate cross
/// product of
///
/// @returns - Returns the vector cross product of u and v
///
//////////////////////////////////////////////////////////////////////////////
        Matrix getCrossProduct(const Matrix &u, const Matrix&v);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function binds segments to closed circuit
///
/// @param Ls - input - list of segments
///
/// @returns - list of polygons
///
//////////////////////////////////////////////////////////////////////////////
        std::list<std::list<Matrix> > cont(std::list<std::list<Matrix> > Ls,double epsilon=1e-2);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief Minimizes polygon
///
/// @param Ls - input - polygon to minimize represented as list of 2D column Vectors
///
/// @returns - Minimized polygon represented as list of 2D column Vectors
///
//////////////////////////////////////////////////////////////////////////////
        std::list<Matrix> min_poly(std::list<Matrix> P,double epsilon=1e-2);


//////////////////////////////////////////////////////////////////////////////
/// @ingroup Debug
///
/// @brief Prints list of matrices
///
/// The only reason of creation is debugin while writing primary functions
///
/// @param fileName - input - name of the file to save list in
/// @param listToPrint - input - self-explanatory
///
/// @returns - Nothing
///
//////////////////////////////////////////////////////////////////////////////
        void printList(const char *fileName,std::list<Matrix> listToPrint);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Debug
///
/// @brief Prints list of lists of matrices
///
/// The only reason of creation is debugin while writing primary functions
///
/// @param fileName - input - name of the file to save list in
/// @param listToPrint - input - self-explanatory
///
/// @returns - Minimized polygon represented as list of 2D column Vectors
///
//////////////////////////////////////////////////////////////////////////////
        void printListOfLists(const char*,std::list<std::list<Matrix> >listToPrint);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function returns matrix of rotation on X
///
/// @param Ls - input - list of segments
///
/// @returns - list of polygons
///
//////////////////////////////////////////////////////////////////////////////
	Matrix rotateX(double angle);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function binds segments to closed circuit
///
/// @param Ls - input - list of segments
///
/// @returns - list of polygons
///
//////////////////////////////////////////////////////////////////////////////
	Matrix rotateY(double angle);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function binds segments to closed circuit
///
/// @param Ls - input - list of segments
///
/// @returns - list of polygons
///
//////////////////////////////////////////////////////////////////////////////
	Matrix rotateZ(double angle);
//////////////////////////////////////////////////////////////////////////////
/// @ingroup Linear algebra
///
/// @brief This function binds segments to closed circuit
///
/// @param Ls - input - list of segments
///
/// @returns - list of polygons
///
//////////////////////////////////////////////////////////////////////////////
	Matrix scale(double x,double y,double z);
	//////////////////////////////////////////////////////////////////////////////
	/// @ingroup Linear algebra
	///
	/// @brief This function binds segments to closed circuit
	///
	/// @param Ls - input - list of segments
	///
	/// @returns - list of polygons
	///
	//////////////////////////////////////////////////////////////////////////////
	Matrix translate(double x,double y,double z);

	//////////////////////////////////////////////////////////////////////////////
	/// @ingroup Linear algebra
	///
	/// @brief This function tests weather point is within a polygon (conv2_test)
	///
	/// @param point - input - Point to test
	/// @param polygon - input - Polygon to test
	///
	/// @returns - 1 if point is within a polygon
	///	-1 if point is outside of a polygon
	///	0 if point is on the border of a polygon
	///
	//////////////////////////////////////////////////////////////////////////////
	int point_is_in_polygon(const Matrix &point, const std::vector<Matrix> &polygon);

	//////////////////////////////////////////////////////////////////////////////
	/// @ingroup Linear algebra
	///
	/// @brief This function tests position of point relating to line
	///
	/// @param A - input - The first point of the line
	/// @param B - input - The second point of the line
	/// @param P - input - A point
	///
	/// @returns - > 0 if point is at the right side of line
	///	< 0 if point is at the left side of line
	///	0 if point is on the line
	///
	//////////////////////////////////////////////////////////////////////////////
	double line_to_point_test(const Matrix &A, const Matrix &B, const Matrix&P);

};


#endif // _Aloschil_Linalg_h
