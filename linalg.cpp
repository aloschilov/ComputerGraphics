///////////////////////////////////////////////////////////////////////////////
/// @file linalg.cpp
/// @brief This file contains implementation of linear algebra modules
///
/// Related Files: The definitions are in linalg.h
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

#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linalg.h"
#include "matrix.h"

using namespace std;

namespace Aloschil
{
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
double ang(const Matrix &v,const Matrix &omega)
{
        if(fabs((v.getTransposed()    *v).getDeterminant())    <EPSILON ||
           fabs((omega.getTransposed()*omega).getDeterminant())<EPSILON)
                return 0.0;
        Matrix d = augment(v,omega);

        double cosValue  = (v.getTransposed()*omega).getDeterminant()
        /(
        sqrt((v.getTransposed()*v).getDeterminant())
        *sqrt((omega.getTransposed()*omega).getDeterminant()));

        if(cosValue < -1) cosValue = -1;
        if(cosValue > 1) cosValue = 1;

        return (d.getDeterminant()>=0 ? 1 : -1)*acos(
        cosValue
        );
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
double atan(const Matrix &p)
{
        return atan2(p.at(1,0),p.at(0,0));
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix
getCrossingParameters(
const Matrix &A,
const Matrix &B,
const Matrix &C,
const Matrix &D)
{
        return !augment(B-A,-(D-C))*(C-A);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix
getCrossPoint(const Matrix &A,const Matrix &B,
                 const Matrix &C,const Matrix &D)
{
        Matrix T = getCrossingParameters(A,B,C,D);
        return A + T[0]*(B-A);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix
getSectionsCrossPoint(const Matrix &A,const Matrix &B,
                 const Matrix &C,const Matrix &D)
{
        Matrix T = getCrossingParameters(A,B,C,D);
        if(!(T[0]>0 && T[0]<1 && T[1]>0 && T[1]<1)) throw NoCrossPoint();
        return A + T[0]*(B-A);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
double
getDistance(const Matrix &A, const Matrix &B)
{
       return sqrt(((B-A).getTransposed()*(B-A)).getDeterminant());
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
SectionsPositionalRelationship
getPositionalRelationship(      const Matrix &A,const Matrix &B,
                                const Matrix &C,const Matrix &D,
                                double rangeOfEquivalence,
                                double rangeOfCollinearity)
{
        if(     fabs(ang(B-A,D-C)/deg)<rangeOfCollinearity ||
                fabs(fabs(ang(B-A,D-C)/deg)-180.0)<rangeOfCollinearity)
        {
                return COLLINEAR;
        }
        else
        {
                Matrix E = getCrossPoint(A,B,C,D);
                char DCBA=0;
                if(getDistance(E,A)<rangeOfEquivalence) DCBA = DCBA | 0x01;
                if(getDistance(E,B)<rangeOfEquivalence) DCBA = DCBA | 0x02;
                if(getDistance(E,C)<rangeOfEquivalence) DCBA = DCBA | 0x04;
                if(getDistance(E,D)<rangeOfEquivalence) DCBA = DCBA | 0x08;

                Matrix T = getCrossingParameters(A,B,C,D);

                if(DCBA==0)
                {// No points are considered to be equal to cross point
                        if(!(T[0]<0 || T[0]>1 || T[1] <0 || T[1]>1))
                                return CROSSING_OR_CROSSED;
                }
		else
                {
                                if(!((DCBA & 0x03) && (DCBA & 0x0C)))
                                { // non-adjacent and non-crossed one
                                        if((DCBA & 0x03) && (T[1]>0 && T[1]<1))
                                        {//touched one
                                                return TOUCHED;
					}
                                        if((DCBA & 0x0C) && (T[0]>0 && T[0]<1))
					{//touching one
                                                return TOUCHING;
					}
                                }
				else
				{// Adjacent by meaning
                                        return ADJACENT;
				}
                }
        }
        return NONINTERSECTING;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
int
random(int max)
{
	int number;
	number = rand();
	return number %  max;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
vector<Matrix>
getRandomPolygon(       const Matrix &T,
                        double a, double b,
                        double teta)
{
        vector<Matrix> resultVector;
        int phi;

        for(phi=0;phi<360;phi+=random(teta))
        resultVector.push_back(T + (a + random(b-a))*Vector(cos(phi*deg),sin(phi*deg)));

        return resultVector;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
vector<Matrix> getRound(
                                const Matrix &a,
                                const Matrix &b,
                                const Matrix &c,
                                double r,
                                double rho,
                                int numberOfSegments,
                                double rangeOfCollinearity)
{
        if(getPositionalRelationship(a,b,b,c,EPSILON,rangeOfCollinearity)==
        COLLINEAR)
        {
                vector<Matrix> result;
                result.push_back(b);
                return result;
        }

        Matrix V = normalize(b-a);
        Matrix W = normalize(c-b);

        double phi = 0.5*ang(V,W);

        double rCandidateAB = rho*getDistance(a,b)/tan(fabs(phi));
        double rCandidateBC = rho*getDistance(b,c)/tan(fabs(phi));

        if(rCandidateAB>=r)
        {
                if(rCandidateBC<r)
                {
                        r = rCandidateBC;
                }

        }
        else
        {
                if(rCandidateAB<rCandidateBC)
                {
                        r = rCandidateAB;
                }
                else
                {
                        r = rCandidateBC;
                }
        }

        Matrix B = sign(phi)* Lrot(V+W);
        Matrix o = b + r/cos(phi)*normalize(B);

        double alfa = atan(sign(phi)*Rrot(V));
        return getArc(o,r,alfa,2*phi,numberOfSegments);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
vector<Matrix> getArc(  const Matrix &center,
                        double radius,
                        double startAngle,
                        double sweepLength,
                        int numberOfSegments)
{
        vector<Matrix> result;

        if(sweepLength>=0)
        {
                if(numberOfSegments<=0) return result;

                double step = sweepLength/numberOfSegments;
                double endAngle = startAngle + sweepLength;

                for(double angle=startAngle;angle<endAngle + EPSILON;angle+=step)
                {
                        double x = radius*cos(angle);
                        double y = radius*sin(angle);
                        Matrix matrixToPush = Vector(x,y)+center;
                        result.push_back(matrixToPush);
                }
        }else
        {
                if(numberOfSegments<=0) return result;

                double step = sweepLength/numberOfSegments;
                double endAngle = startAngle + sweepLength;

                for(double angle=startAngle;angle>=endAngle - EPSILON;angle+=step)
                {
                        double x = radius*cos(angle);
                        double y = radius*sin(angle);
                        Matrix matrixToPush = Vector(x,y)+center;
                        result.push_back(matrixToPush);
                }

        }

        return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix Lrot(const Matrix &vectorToRotate)
{
        Matrix result;
        result.setSize(2,1);
        result[0] = -vectorToRotate.at(1,0);
        result[1] = vectorToRotate.at(0,0);
        return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix Rrot(const Matrix &vectorToRotate)
{
        Matrix result;
        result.setSize(2,1);
        result[0] = vectorToRotate.at(1,0);
        result[1] = -vectorToRotate.at(0,0);
        return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix normalize(const Matrix &vectorToNormalize)
{
        double length =sqrt( (vectorToNormalize.getTransposed()*vectorToNormalize).
                getDeterminant());
        return (1/length)*vectorToNormalize;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
double sign(double x)
{
        if(x==0) return 0;
        if(x>0) return 1;
        return -1;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
std::vector<std::vector<Matrix> > getClosedCircuits(
                                std::vector<std::vector<Matrix> > listOfSegments,
                                double rangeOfEquivalence)
{
        size_t numberOfSegments = listOfSegments.size();
        vector<vector<Matrix> > emptyVectorOfPolygons;
        vector<vector<Matrix> > result;
        vector<Matrix> ab;
        vector<Matrix > curCircuit;
        size_t currentCircuitSize = 0;

        bool newCircuitJustStarted = true;
        bool lastCircuitJustClosed = true;

        if(listOfSegments.size()<3)
                return emptyVectorOfPolygons;

        vector<bool> notUsed(listOfSegments.size(),true);

        bool thereAreChanges = true;

        while(true)
        {
                if(!thereAreChanges) return emptyVectorOfPolygons;

                thereAreChanges = false;

                for(size_t i=0;i<numberOfSegments;++i)
                {
                        if(find(notUsed.begin(),notUsed.end(),true)==notUsed.end())
                        {
                                if(lastCircuitJustClosed)
                                        return result;
                                else
                                        return emptyVectorOfPolygons;
                        }
                        lastCircuitJustClosed = false;
                        if(notUsed[i])
                        {
                                ab = listOfSegments[i];
                                if(newCircuitJustStarted)
                                {
                                        curCircuit.clear();
                                        curCircuit.push_back(ab[0]);
                                        curCircuit.push_back(ab[1]);
                                        currentCircuitSize = 2;
                                        newCircuitJustStarted = false;
                                        notUsed[i] = false;
                                        thereAreChanges = true;
                                }
                                else
                                {
                                        if(getDistance(
                                        curCircuit[currentCircuitSize-1],
                                        ab[0])<rangeOfEquivalence)
                                        {
                                                curCircuit.push_back(ab[1]);
                                                ++currentCircuitSize;
                                                notUsed[i] = false;
                                                thereAreChanges = true;
                                        }
                                        else
                                        {
                                               if(getDistance(
                                                curCircuit[currentCircuitSize-1],
                                                ab[1])<rangeOfEquivalence)
                                                {
                                                        curCircuit.push_back(ab[0]);
                                                        ++currentCircuitSize;
                                                        notUsed[i] = false;
                                                        thereAreChanges = true;
                                                }
                                                else
                                                {
                                                        continue;
                                                }
                                        }

                                        if(getDistance(curCircuit[currentCircuitSize-1],
                                        curCircuit[0])<rangeOfEquivalence && currentCircuitSize>3)
                                        {
                                                result.push_back(curCircuit);
                                                newCircuitJustStarted = true;
                                                lastCircuitJustClosed = true;
                                        }
                                        else
                                        {
                                                continue;
                                        }

                                }
                        }
                        else
                        {
                                continue;
                        }
                }
        }
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix getCrossProduct(const Matrix &u, const Matrix&v)
{
        Matrix result;
        result.setSize(3,1);
        result[0]=u.at(1,0)*v.at(2,0) - v.at(1,0)*u.at(2,0);
        result[1]=-(u.at(0,0)*v.at(2,0) - v.at(0,0)*u.at(2,0));
        result[2]=u.at(0,0)*v.at(1,0) - v.at(0,0)*u.at(1,0);
        return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
list<list<Matrix> > cont(list<list<Matrix> > Ls,double epsilon)
{
        list<list<Matrix> > LP;
        list<list<Matrix> > S;
        list<Matrix> ab;
        list<Matrix> C;

        printListOfLists("Ls.list",Ls);

        size_Ls:
        if(!Ls.empty())
        {
                S = Ls;
                printListOfLists("S.list",S);
                ab_Ls:
                ab = (*Ls.begin());
                printList("ab.list",ab);
                if(C.empty())
                {
                        C=ab;
                        printList("C.list",C);
                        goto Ls_Ls_ab;
                }
                else
                {
                        if(getDistance(C.back(),ab.front())<=epsilon)
                        {
                                C.push_back(ab.back());
                                goto Ls_Ls_ab;
                        }
                        else
                        {
                                if(getDistance(C.back(),ab.back())<=epsilon)
                                {
                                        C.push_back(ab.front());
                                        goto Ls_Ls_ab;
                                }
                                else
                                {
                                        list<Matrix> tmp = Ls.front();
                                        Ls.erase(Ls.begin());
                                        Ls.push_back(tmp);
                                        printListOfLists("Ls.list",Ls);
                                        // ���� �������� �������������� ����, ��� ����� ������� {Ls=S}
                                        if(equal(S.begin(),S.end(),Ls.begin()))
                                        {
                                                        list<Matrix>::iterator endSegment = C.end();
                                                        --endSegment;
                                                        C.erase(endSegment);
                                                        goto size_Ls;
                                        }
                                        else
                                        {
                                                goto ab_Ls;
                                        }
                                }
                        }
                }

                Ls_Ls_ab:
                // �������� ������� �� ������
                list<list<Matrix> >::iterator curSegment = Ls.begin();
                list<list<Matrix> >::iterator endSegment = Ls.end();
                // ����� ���� �������� � ������
                for(;curSegment!=endSegment;++curSegment)
                {
                        if(     (*curSegment).front()==ab.front() &&
                                (*curSegment).back()==ab.back())
                        {
                                Ls.erase(curSegment);
                                break;
                        }
                }
                printListOfLists("Ls.list",Ls);
                if(C.size()<4) goto size_Ls;


                printList("C.list",C);
                list<Matrix>::iterator curVertex = C.begin();
                list<Matrix>::iterator endVertex = C.end();
                list<Matrix>::iterator C_m = C.end();--C_m;
                --endVertex;
                --endVertex;
                --endVertex;

                for(;curVertex!=endVertex;++curVertex)
                {
                        if(getDistance(*curVertex,*C_m)>epsilon)
                                continue;
                        if(curVertex!=C.begin())
                        {
                                list<Matrix>::iterator k = C.begin();
                                for(;k!=curVertex;++k)
                                {
                                        list<Matrix > c1c2;
                                        c1c2.push_back(C.front());
                                        C.erase(C.begin());
                                        c1c2.push_back(C.front());
                                        Ls.push_back(c1c2);
                                }
                        }
                        printList("C.list",C);
                        LP.push_back(min_poly(C));
                        printListOfLists("LP.list",LP);
                        printListOfLists("Ls.list",Ls);
                        goto size_Ls;

                }
                goto size_Ls;
        }
        printListOfLists("LP.list",LP);
        return LP;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
list<Matrix> min_poly(list<Matrix> P,double epsilon)
{
        list<Matrix>::iterator pm = P.begin();
        list<Matrix>::iterator i = P.begin();
        list<Matrix>::iterator pEnd = P.end();
        double t;
        for(++i;i!=pEnd;++i)
        {
                Matrix V = (*i) - (*pm);
                //double abs = sqrt((V.getTransposed()*V).getDeterminant());
                if(sqrt((V.getTransposed()*V).getDeterminant())>epsilon)
                {
                        if(pm!=P.begin())
                        {
                                list<Matrix>::iterator pmM1 = pm;
                                --pmM1;
                                Matrix W = (*pm) - (*pmM1);
                                t = (W.getTransposed()*(V + W)).getDeterminant() / (W.getTransposed()*W).getDeterminant();
                                if(fabs(((normalize(W).getTransposed())*normalize(V)).getDeterminant())>=1)
                                {
                                        if(t>1)
                                        {//(1,inf)
                                                (*pm) = (*i);
                                        }
                                        else
                                        {
                                                if(t<=0)
                                                {// (-inf,0]
                                                        (*pmM1)=(*pm);
                                                        (*pm) = (*i);
                                                }
                                                else
                                                {// (0,1]
                                                        continue;
                                                }
                                        }
                                }
                                else
                                {
                                        ++pm;
                                        (*pm) = (*i);
                                }
                        }
                        else
                        {
                                ++pm;
                                (*pm) = (*i);
                        }
                }
        }

        list<Matrix>::iterator pmM1=pm;
        --pmM1;
        list<Matrix>::iterator pmP1=pm;
        ++pmP1;

        list<Matrix>::iterator p2 = P.begin();
        ++p2;

        if(fabs((normalize(P.front() - (*pmM1))*normalize((*p2) - P.front()).getTransposed()).getDeterminant())<1)
        {
                P.erase(pmP1,P.end());
        }
        else
        {
                P.erase(P.begin());
                P.erase(pm,P.end());
                P.push_back(P.front());
        }
        return P;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
void printList(const char *fileName,list<Matrix> listToPrint)
{
        FILE *file;
        if((file=fopen(fileName,"w"))==NULL)
        {
         return;
        }
        list<Matrix>::iterator curMatrix = listToPrint.begin();
        list<Matrix>::iterator endMatrix = listToPrint.end();

        for(;curMatrix!=endMatrix;++curMatrix)
        {
                for(unsigned int i=0;i<(*curMatrix).y;i++)
                {
                        for(unsigned int j=0;j<(*curMatrix).x;j++)
                        {
                                fprintf(file,"|%10.3g",
                                (*curMatrix).field[i*(*curMatrix).x+j]);
                        }
                        fprintf(file,"\n");
                }
                fprintf(file,"\n");
        }
        fclose(file);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
void printListOfLists(const char *fileName,list<list<Matrix> >listToPrint)
{
        FILE *file;
        if((file=fopen(fileName,"w"))==NULL)
        {
         return;
        }
        list<list<Matrix> >::iterator curList = listToPrint.begin();
        list<list<Matrix> >::iterator endList = listToPrint.end();

        for(;curList!=endList;++curList)
        {
                list<Matrix>::iterator curMatrix = (*curList).begin();
                list<Matrix>::iterator endMatrix = (*curList).end();

                for(;curMatrix!=endMatrix;++curMatrix)
                {
                        for(unsigned int i=0;i<(*curMatrix).y;i++)
                        {
                                for(unsigned int j=0;j<(*curMatrix).x;j++)
                                {
                                        fprintf(file,"|%10.3g",
                                        (*curMatrix).field[i*(*curMatrix).x+j]);
                                }
                                fprintf(file,"\n");
                        }
                        fprintf(file,"\n");
                }
        }
        fclose(file);
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix rotateX(double angle)
{
	Matrix tmp;
	tmp.setSize(4, 4);
	tmp.at(0,0) = 1.0;
	tmp.at(1,1) = cos(angle);
	tmp.at(1,2) = sin(angle);
	tmp.at(2,1) = -sin(angle);
	tmp.at(2,2) = cos(angle);
	tmp.at(3,3) = 1.0;
	return tmp.getTransposed();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix rotateY(double angle)
{
	Matrix tmp;
	tmp.setSize(4, 4);
	tmp.at(0,0) = cos(angle);
	tmp.at(0,2) = -sin(angle);
	tmp.at(1,1) = 1.0;
	tmp.at(2,0) = sin(angle);
	tmp.at(2,2) = cos(angle);
	tmp.at(3,3) = 1.0;
	return tmp.getTransposed();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix rotateZ(double angle)
{
	Matrix tmp;
	tmp.setSize(4, 4);
	tmp.at(0,0) = cos(angle);
	tmp.at(0,1) = sin(angle);
	tmp.at(1,0) = -sin(angle);
	tmp.at(1,1) = cos(angle);
	tmp.at(2,2) = 1.0;
	tmp.at(3,3) = 1.0;
	return tmp.getTransposed();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix scale(double x,double y,double z)
{
	Matrix tmp;
	tmp.setSize(4,4);
	tmp.at(0,0)=x;
	tmp.at(1,1)=y;
	tmp.at(2,2)=z;
	tmp.at(3,3)=1.0;
	return tmp.getTransposed();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
Matrix translate(double x,double y,double z)
{
	Matrix tmp;
	tmp.setSize(4,4);
	tmp.at(0,0)=1.0;
	tmp.at(1,1)=1.0;
	tmp.at(2,2)=1.0;
	tmp.at(3,3)=1.0;
	tmp.at(3,0)=x;
	tmp.at(3,1)=y;
	tmp.at(3,2)=z;
	return tmp.getTransposed();
}
double _M_PI = 3.141592654;
double const deg = _M_PI/180.0;

//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
int point_is_in_polygon(const Matrix &point, const std::vector<Matrix> &polygon) {
	int r = 0;
	int l = 0;
	int e = 0;
	size_t segments_count = polygon.size();

	for(size_t current_segment=0; current_segment<segments_count; ++current_segment) {
		double point_test_result = line_to_point_test(
				polygon[current_segment],
				polygon[(current_segment+1)%segments_count],
				point);
		if(point_test_result<0) {
			l = 1;
		} else if (point_test_result>0) {
			r = 1;
		} else {
			e = 1;
		}

		if(l*r != 0) {
			return -1;
		}
	}
	return -(e-1);
}

//----------------------------------------------------------------------------
// This is a public API. Refer to linalg.h for details.
//----------------------------------------------------------------------------
double line_to_point_test(const Matrix &A, const Matrix &B, const Matrix&P) {
	Matrix M;
	M.setSize(2, 2);
	M.at(0, 0) = (P - A)[0];
	M.at(1, 0) = (P - A)[1];

	M.at(0, 1) = (B - A)[0];
	M.at(1, 1) = (B - A)[1];

	return M.getDeterminant();
}

} // namespace Aloschil
