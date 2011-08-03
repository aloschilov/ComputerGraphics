#include "linalg.h"

#include "object3D.h"

namespace Aloschil
{
using namespace std;
//----------------------------------------------------------------------------
// This is a public API. Refer to object3D.h for details.
//----------------------------------------------------------------------------
void Face::setPoints(std::vector<Matrix> pointsToSet)
{
    points=pointsToSet;
	updateNormal();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to object3D.h for details.
//----------------------------------------------------------------------------
Face Face::getInvertedFace()
{
	Face result;
	int numberOfPoints = (int)points.size();
	for(int i=numberOfPoints-1;i>=0;--i)
	{
		result.points.push_back(points[i]);
	}
	result.updateNormal();
	return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to object3D.h for details.
//----------------------------------------------------------------------------
void Face::print()
{
	vector<Matrix>::iterator curPoint=points.begin();
	vector<Matrix>::iterator endPoint=points.end();
	printf("Face points:\n");
	for(;curPoint!=endPoint;++curPoint)
	{
		(*curPoint).print();
	}
	printf("Face normal:\n");
	normal.print();
}
//----------------------------------------------------------------------------
// This is a public API. Refer to object3D.h for details.
//----------------------------------------------------------------------------
Face operator *(const Matrix &m,const Face &face)
{
	Face result;
	vector<Matrix>::const_iterator curPoint=face.points.begin();
	vector<Matrix>::const_iterator endPoint=face.points.end();
	for(;curPoint!=endPoint;++curPoint)
	{
		result.points.push_back(m*(*curPoint));
	}
	result.updateNormal();
	return result;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to object3D.h for details.
//----------------------------------------------------------------------------
void Face::updateNormal()
{
	normal = normalize(getCrossProduct(points[1]-points[0],points[2]-points[1]));
}
}
