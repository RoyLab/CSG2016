// (c) 2010 Geom e.U. Bernhard Kornberger, Graz/Austria. All rights reserved.
//
// This file is part of the Fade2D library. You can use it for your personal
// non-commercial research. Licensees holding a commercial license may use this 
// file in accordance with the Commercial License Agreement provided 
// with the Software.
//
// This software is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING 
// THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are not clear 
// to you.
// 
// Author: Bernhard Kornberger, bkorn (at) geom.at 
// http://www.geom.at


#pragma once
#if GEOM_PSEUDO3D==GEOM_TRUE

#include "Segment2.h"
#include <map>

#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE 
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined 
#endif 


struct TriNode; // FWD
class Triangle2; // FWD
class GeomTest; // FWD
class Segment2; // FWD

typedef std::multimap<Point2,Segment2*> MMPS;
typedef MMPS::iterator MMPSIt;

class IsoContours  
{

public:
/**  
*
* Constructs a fast datastructure suited to compute intersections of a 
* horizontal plane with the given triangles. Typically only one such object
* is made and then its getContours() function is called several times.
*   
*/

	CLASS_DECLSPEC
	IsoContours(std::vector<Triangle2*>& vTriangles);
	
	CLASS_DECLSPEC
	~IsoContours();
/** 
 * 
 * Returns the smallest z coordinate
 */
	CLASS_DECLSPEC
	double getMinHeight();
/**  
 *
 * Returns the largest z-coordinate
 */
	CLASS_DECLSPEC
	double getMaxHeight();
/**   
 *
 * Computes the intersection of a horizontal plane at a certain height with all triangles and
 * returns a vector of assembled polygons and polylines.
 */
	CLASS_DECLSPEC
	bool getContours(double height,std::vector<std::vector<Segment2> >& vvContours,bool bVerbose);
protected:
	CLASS_DECLSPEC
	bool getIntersectionSegments(double height,std::vector<Segment2*>& vIntersectionSegments);
	void createContours(MMPS& mmPS,std::vector<std::vector<Segment2> >& vvContours);
	void getIntersection(Triangle2* pT,double height,int i,std::vector<Point2>& vIntersectionPoints);

	TriNode* pTree;
	std::vector<Triangle2*> vTriangles;
	GeomTest* pGeomPredicates;
	std::set<double> sForbiddenHeights;
};


} // NAMESPACE FADE2D


#endif


