#include <stdlib.h>

#include "linalg.h"
#include "fractals.h"

namespace Aloschil
{
using namespace std;
//----------------------------------------------------------------------------
// This is a public API. Refer to fractals.h for details.
//----------------------------------------------------------------------------
SimpleGeometricFractal::SimpleGeometricFractal()
{
	deletionType = DeletePrevious;
	overallPolicy= DontMirrorOverallPolicy;
	onIterationPolicy= DontMirrorOnIterationPolicy;
	onSegmentNumberPolicy= DontMirrorOnSegmentNumberPolicy;
}
//----------------------------------------------------------------------------
// This is a public API. Refer to fractals.h for details.
//----------------------------------------------------------------------------
void SimpleGeometricFractal::calculate(const list<Segment> &_iterator,const list<Segment> &_initializer,unsigned int numberOfIterations,double delta)
{
	iterator=_iterator;
	
	iterations.clear();

	iterations.push_back(_initializer);

	double overallYScaleSign = (overallPolicy==DontMirrorOverallPolicy) ? 1 : -1;

	for(unsigned int i=0;i<numberOfIterations;++i)
	{
		list<Segment> currIterationSegments;
		list<Segment>::iterator curSegment = iterations[i].begin();
		list<Segment>::iterator endSegment = iterations[i].end();

		double onIterationScaleSign;

		switch(onIterationPolicy)
		{
			case MirrorOnEvenOnIterationPolicy:
				onIterationScaleSign = (i%2==0) ? -1 : 1;
				break;
			case MirrorOnOddOnIterationPolicy:
				onIterationScaleSign = (i%2!=0) ? -1 : 1;
				break;
			case DontMirrorOnIterationPolicy:
				onIterationScaleSign = 1;
				break;
		}

		unsigned int segmentNumber=0;

		for(;curSegment!=endSegment;++curSegment,++segmentNumber)
		{
			double onSegmentNumberScaleSign;
			switch(onSegmentNumberPolicy)
			{
				case MirrorOnEvenOnSegmentNumberPolicy:
					onSegmentNumberScaleSign = (segmentNumber%2==0) ? -1 : 1;
					break;
				case MirrorOnOddOnSegmentNumberPolicy:
					onSegmentNumberScaleSign = (segmentNumber%2!=0) ? -1 : 1;
					break;
				case DontMirrorOnSegmentNumberPolicy:
					onSegmentNumberScaleSign = 1;
					break;
			}

			double curSign = onSegmentNumberScaleSign*onIterationScaleSign*overallYScaleSign;
			
			double curSegmentLength = (*curSegment).getLength();
			if(curSegmentLength>delta)
			{
				Matrix curSegmentVector = (*curSegment).getVector();
				Matrix &curFirstVertex = (*curSegment).p1;

				Matrix transformation = translate(curFirstVertex[0], curFirstVertex[1], curFirstVertex[2])*
				rotateZ(ang(Vector(1.0,0.0),Vector(curSegmentVector[0],curSegmentVector[1])))*
				scale(curSegmentLength,curSign*curSegmentLength,curSegmentLength);

				list<Segment>::iterator curIteratorSegment = iterator.begin();
				list<Segment>::iterator endIteratorSegment = iterator.end();

				for(;curIteratorSegment!=endIteratorSegment ;++curIteratorSegment)
				{
					Matrix tmp1 = transformation*(*curIteratorSegment).p1;
					Matrix tmp2 = transformation*(*curIteratorSegment).p2;
					Segment tmp;
					tmp.p1 = tmp1;
					tmp.p2 = tmp2;
					currIterationSegments.push_back(tmp);
				}

				switch(deletionType)
				{
					case KeepPrevious:
						currIterationSegments.push_back(*curSegment);
						break;
					case DeletePrevious:
						break;
					case DeleteOnEvenIterations:
						if(i%2)
							currIterationSegments.push_back(*curSegment);						
						break;
					case RandomDeletion:
						if(rand()>RAND_MAX/2.0)
							currIterationSegments.push_back(*curSegment);
						break;
				}
			}
			else
			{
				currIterationSegments.push_back(*curSegment);
			}
		}
		iterations.push_back(currIterationSegments);
	}
}

void IterableFunctionsSystemFractal::calculate(
	const vector<Matrix> &baseVertex,
	const vector<Matrix> &transformationTriples,
	unsigned int numberOfIterations)
{
	iterations.clear();
	// Calculating transformation matrices

	vector<Matrix> transformationMatrices;

	if(baseVertex.empty()) return;

	Matrix baseTrianlgeMatrix = baseVertex[0];
	baseTrianlgeMatrix = augment(baseTrianlgeMatrix, baseVertex[1]);
	baseTrianlgeMatrix = augment(baseTrianlgeMatrix, baseVertex[2]);

	Matrix inverseBaseTriangleMatrix = (!baseTrianlgeMatrix);

	unsigned int numberOfTransformations = transformationTriples.size()/3;
	for(unsigned int i=0;i<numberOfTransformations;++i)
	{
		Matrix endTrianlgeMatrix = transformationTriples[i*3];
		endTrianlgeMatrix = augment(endTrianlgeMatrix, transformationTriples[i*3+1]);
		endTrianlgeMatrix = augment(endTrianlgeMatrix, transformationTriples[i*3+2]);
		transformationMatrices.push_back(endTrianlgeMatrix*inverseBaseTriangleMatrix);
	}

	iterations.push_back(transformationTriples);
	
	for(unsigned int i=1;i<numberOfIterations;++i)
	{
		vector<Matrix> currentIterationVertices;
		unsigned int numberOfVerticesInPrIt = iterations[i-1].size();
		for(unsigned int curVertex=0;curVertex<numberOfVerticesInPrIt;++curVertex)
		{
			unsigned int numberOfTransformations = transformationMatrices.size();
			for(unsigned int curTransform=0;curTransform<numberOfTransformations;++curTransform)
				currentIterationVertices.push_back(transformationMatrices[curTransform]*iterations[i-1][curVertex]);
		}
		iterations.push_back(currentIterationVertices);
	}
}
}
