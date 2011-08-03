#ifndef _Aloschil_Fractals_h_
#define _Aloschil_Fractals_h_

///////////////////////////////////////////////////////////////////////////////
/// @file fractals.h
/// @brief This file contains declaration of fractal series classes
///
/// Related Files: The implementation is in fractals.cpp
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
#include <list>
#include <vector>

#include "Matrix.h"

namespace Aloschil
{
struct Segment
{
	Matrix p1;
	Matrix p2;
	double getLength()
	{
		double p1_x=p1[0]/p1[3];
		double p1_y=p1[1]/p1[3];
		double p1_z=p1[2]/p1[3];

		double p2_x=p2[0]/p2[3];
		double p2_y=p2[1]/p2[3];
		double p2_z=p2[2]/p2[3];

		double deltaX = p2_x-p1_x;
		double deltaY = p2_y-p1_y;
		double deltaZ = p2_z-p1_z;

		return sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
	}
	Matrix getVector()
	{
		double p1_x=p1[0]/p1[3];
		double p1_y=p1[1]/p1[3];
		double p1_z=p1[2]/p1[3];

		double p2_x=p2[0]/p2[3];
		double p2_y=p2[1]/p2[3];
		double p2_z=p2[2]/p2[3];

		double deltaX = p2_x-p1_x;
		double deltaY = p2_y-p1_y;
		double deltaZ = p2_z-p1_z;
		return Vector4(deltaX,deltaY,deltaZ);
	}
};

class SimpleGeometricFractal
{
	public:
		enum DeletionType
		{
			KeepPrevious,
			DeletePrevious,
			DeleteOnEvenIterations,
			RandomDeletion
		};
		enum OverallPolicy
		{
			MirrorOverallPolicy,
			DontMirrorOverallPolicy
		};
		enum OnIterationPolicy
		{
			MirrorOnEvenOnIterationPolicy,
			MirrorOnOddOnIterationPolicy,
			DontMirrorOnIterationPolicy
		};
		enum OnSegmentNumberPolicy
		{
			MirrorOnEvenOnSegmentNumberPolicy,
			MirrorOnOddOnSegmentNumberPolicy,
			DontMirrorOnSegmentNumberPolicy
		};
		
		SimpleGeometricFractal();
		void calculate(const std::list<Segment> &_iterator,const std::list<Segment> &_initializer,unsigned int numberOfIterations,double delta=0.5);
		std::list<Segment> & operator [] (unsigned int i) {return iterations[i];}
		std::list<Segment> & at (unsigned int i) {return iterations.at(i);}

		void setDeletionType(DeletionType newDeletionType){deletionType=newDeletionType;}		
		void setOverallPolicy(OverallPolicy newOverallPolicy){overallPolicy=newOverallPolicy;}
		void setOnIterationPolicy(OnIterationPolicy newOnIterationPolicy){onIterationPolicy=newOnIterationPolicy;}
		void setOnSegmentNumberPolicy(OnSegmentNumberPolicy newOnSegmentNumberPolicy){onSegmentNumberPolicy=newOnSegmentNumberPolicy;}
	private:
		std::list<Segment> iterator;
		std::list<Segment> initializer;
		
		std::vector<std::list<Segment> > iterations;
		
		DeletionType deletionType;
		OverallPolicy overallPolicy;
		OnIterationPolicy onIterationPolicy;
		OnSegmentNumberPolicy onSegmentNumberPolicy;

};

class IterableFunctionsSystemFractal
{
	public:
		IterableFunctionsSystemFractal(){}
		void calculate(const std::vector<Matrix> &baseVertex,const std::vector<Matrix> &transformationTriples,unsigned int numberOfIterations);
		std::vector<std::vector<Matrix> > & getIterations() {return iterations;}
	private:
		std::vector<std::vector<Matrix> > iterations;
};
}
#endif // _Aloschil_Fractals_h_

