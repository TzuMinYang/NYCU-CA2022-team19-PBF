#include <MathLib/Segment.h>
#include <GLFW/glfw3.h>

/**
	Constructor
*/
Segment::Segment(const P3D& a_, const P3D& b_){
	this->a = a_;
	this->b = b_;
}

/**
	Destructor
*/
Segment::~Segment(){
	
}
	
/**
	Copy constructor
*/
Segment::Segment(const Segment& other){
	this->a = other.a;
	this->b = other.b;
}

/**
	Copy operator
*/
Segment& Segment::operator = (const Segment& other){
	this->a = other.a;
	this->b = other.b;

	return *this;
}

/**
	Draws the segment - for visualization purposes.
*/
void Segment::draw(){
    glPointSize(4);
	glBegin(GL_POINTS);
		glVertex3d(a[0], a[1], a[2]);
		glVertex3d(b[0], b[1], b[2]);
	glEnd();

	glBegin(GL_LINES);
		glVertex3d(a[0], a[1], a[2]);
		glVertex3d(b[0], b[1], b[2]);
	glEnd();

	glPointSize(1);
}

/**
	This method returns the point on the current segment that is closest to the point c that is passed in as a paremter.
*/
P3D Segment::getClosestPointTo(const P3D& c) const{
	//we'll return the point d that belongs to the segment, such that: cd . ab = 0
	//BUT, we have to make sure that this point belongs to the segment. Otherwise, we'll return the segment's end point that is closest to D.
	V3D v(a, b);

	double len_squared = v.dot(v);
	//if a==b, it means either of the points can qualify as the closest point
	if (IS_ZERO(len_squared)){
		return a;
	}
	
	double mu = V3D(a, c).dot(v) / len_squared;
	boundToRange(&mu, 0, 1);
	//the point d is at: a + mu * ab
	return a + v * mu;
}

/**
	This method returns the segment that connects the closest pair of points - one on the current segment, and one on the segment that is passed in. The
	'a' point of the resulting segment will lie on the current segment, while the 'b' point lies on the segment that is passed in as a parameter.
*/
Segment Segment::getShortestSegmentTo(const Segment& other) const{
	//MAIN IDEA: the resulting segment should be perpendicular to both of the original segments. Let point c belong to the current segment, and d belong to
	//the other segment. Then a1b1.cd = 0 and a2b2.cd = 0. From these two equations with two unknowns, we need to get the mu_c and mu_d parameters
	//that will let us compute the points c and d. Of course, we need to make sure that they lie on the segments, or adjust the result if they don't.

	//unfortunately, there are quite a few cases we need to take care of. Here it is:
	V3D tmp1(this->a, other.a), tmp2(this->a, this->b), tmp3(other.a, other.b);

	double A = tmp1.dot(tmp2);
	double B = tmp2.dot(tmp3);
	double C = tmp2.dot(tmp2);
	double D = tmp3.dot(tmp3);
	double E = tmp1.dot(tmp3);

	//now a few special cases:
	if (IS_ZERO(C)){
		//current segment has 0 length
		return Segment(a, other.getClosestPointTo(this->a));
	}
	if (IS_ZERO(D)){
		//other segment has 0 length
		return Segment(getClosestPointTo(other.a), other.a);
	}

	if (IS_ZERO(C*D - B*B)){
		//this means that the two segments are coplanar and parallel. In this case, there are
		//multiple segments that are perpendicular to the two segments (lines before the truncation really).

		//we need to get the projection of the other segment's end point on the current segment
		double mu_a2 = tmp1.dot(tmp2) / (tmp2.dot(tmp2));
		double mu_b2 = V3D(a, other.b).dot(tmp2) / (tmp2.dot(tmp2));

		//we are now interested in the parts of the segments that are in common between the two input segments
		boundToRange(&mu_a2, 0, 1);
		boundToRange(&mu_b2, 0, 1);

		//closest point on the current segment will lie at the midpoint of mu_a2 and mu_b2
		return Segment(this->a + tmp2 * ((mu_a2 + mu_b2) / 2.0), other.getClosestPointTo(this->a + tmp2 * ((mu_a2 + mu_b2) / 2.0)));
	}

	//ok, now we'll find the general solution for two lines in space:
	double mu_c = (A*D - E*B) / (C*D-B*B);
	double mu_d = (mu_c*B - E) / D;

	//if the D point or the C point lie outside their respective segments, clamp the values
	boundToRange(&mu_c, 0, 1);
	boundToRange(&mu_d, 0, 1);

	return Segment(this->a + tmp2 * mu_c, other.a + tmp3 * mu_d);
}

