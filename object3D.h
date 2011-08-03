#ifndef _object3D_Aloschil_h_
#define _object3D_Aloschil_h_

#include <vector>
#include "matrix.h"

namespace Aloschil
{
class Face
{
		std::vector<Matrix> points;
		Matrix normal;
	public:
		Face(){}
//////////////////////////////////////////////////////////////////////////////
/// @brief This function sets points of the face
///
/// Points are set in order of traversal and normal is also calculated
/// normal is calculated by first three points.
/// Default traversal is counter-clockwise
///
/// @returns – nothing
///
//////////////////////////////////////////////////////////////////////////////
		void setPoints(std::vector<Matrix> pointsToSet);
//////////////////////////////////////////////////////////////////////////////
/// @brief This function returns face with the same points set but inverted normal
///
/// Note: traversal of face's points is also inverted. <A HREF="#faceInversionImage">See image</A>
///  <A NAME="faceInversion"> 
///  <IMG src=images/faceInversion.jpg ALIGN=left>
///  </A>
///
/// @returns – Return codes and their meanings.
///
//////////////////////////////////////////////////////////////////////////////
		Face getInvertedFace();

		std::vector<Matrix> &getPoints(){return points;}
		Matrix & getNormal(){return normal;}
		Matrix & operator [](int i){return points[i];}


		friend Face operator *(const Matrix &m,const Face &face);
		void print();
		void updateNormal();
};

class Object3D
{
		std::vector<Face> faces;
	public:
		Object3D(){}
		std::vector<Face> & getFaces(){return faces;}
		~Object3D(){}

};
}
#endif /*_object3D_Aloschil_h_*/
