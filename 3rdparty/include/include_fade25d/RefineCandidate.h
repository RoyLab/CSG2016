#pragma once


#include "Point2.h"
#include "Fade_2D.h"
#include "RefineParams.h"

#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

class Triangle2; // Fwd
class Dt2; // Fwd
class ConstraintSegment2; // Fwd

class RefineCandidate
{
public:
	RefineCandidate(Dt2* pDt_,Triangle2* pTriangle_);

	Triangle2* getTriangle() const;
	Point2 getCircumcenter() const;
	double getRadius() const;
	std::pair<Triangle2*,int> getContainingTriangleAndIndex() const;
	Triangle2* getContainingTriangle() const;
	bool dualCloseToInfinity() const;
	bool isConstraint(int ith) const;
	ConstraintSegment2* getConstraintSegment(int ith) const;
	bool operator>(const RefineCandidate& pOther) const;
	bool isAlive() const;
	int debug_label;
	static int debug_runningLabel;
	void debug();
	int getSplitIndex() const;
	void computeMinMax();
	double getQuality2D() const;
	double getMinSqSideLength2D() const;
	double getMaxSqSideLength2D() const;
	double getAspect2D();
#if GEOM_PSEUDO3D==GEOM_TRUE
	double getQuality25D() const;
	double getMinSqSideLength25D() const;
	double getMaxSqSideLength25D() const;
	double getAspect25D();
#endif

	bool hasNarrowConstraintAngle(double requiredAngleDegree) const;
	bool needRefinementDueToGrowth(double growLimit,double growFactorMinArea );

#if GEOM_PSEUDO3D==GEOM_TRUE
	bool needRefinementDueToSlope(double slopeLimit,Fade_2D* pHeightGuideTriangulation) const;
	bool needRefinementDueToHeightError(Fade_2D* pHeightGuideTriangulation,double maxHeightError) const;
#endif

	bool needRefinement( RefineParams& refineParams );


protected:

	int walkUntilConstraint(Point2* pVtx,Triangle2*& pStartTriangle);
	Dt2* pDt;
	Triangle2* pTriangle;
	Triangle2* pContainingT;
	Point2 circumcenter;
	double radius;
	int splitIndex;
	bool bAccurateCC;
	Point2* vOriginalPoints[6];
	std::vector<double> vSqLengthsIgnoringZ;

#if GEOM_PSEUDO3D==GEOM_TRUE
	std::vector<double> vSqLengthsConsideringZ;
#endif


};


} // Namespace

