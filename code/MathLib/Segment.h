#pragma once

#include <MathLib/P3D.h>
#include <MathLib/V3D.h>

/*==========================================================================================================================================*
 *	This class provides the implementation of a segment, and a few useful method associated with them.                                      *
 *==========================================================================================================================================*/
class Segment{
public:
	//the end points of the segment
	P3D a, b;

	Segment(const P3D& a_, const P3D& b_);
	
	Segment(){
		a = P3D();
		b = P3D();
	}
	
	/**
		Copy constructor
	*/
	Segment(const Segment& other);

	/**
		Copy operator
	*/
	Segment& operator = (const Segment& other);


	~Segment();

	/**
		This method returns the point on the current segment that is closest to the point passed in as a paremter.
	*/
	P3D getClosestPointTo(const P3D& c) const;

	/**
		This method returns the segment that connects the closest pair of points - one on the current segment, and one on the segment that is passed in. The
		'a' point of the resulting segment will lie on the current segment, while the 'b' point lies on the segment that is passed in as a parameter.
	*/
	Segment getShortestSegmentTo(const Segment& other) const;

	/**
		Draws the segment - for visualization purposes.
	*/
	void draw();

};



