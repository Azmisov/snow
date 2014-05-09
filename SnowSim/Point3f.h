/**
 * A 3-dimensional data value
 */

#ifndef POINT3F_H
#define	POINT3F_H

class Point3f {
public:
	//Variables
	float loc[3];
	
	//Constructors
	Point3f();
	Point3f(float val);
	Point3f(float x, float y, float z);
	Point3f(const Point3f& orig);
	virtual ~Point3f();
	
	void normalize();
	float dot(const Point3f &v) const;
	float sum() const;
	float length() const;
	
	//Operator Overloads
	//Vector * Scalar
	const Point3f operator*(const float& c) const;
	Point3f& operator*=(const float& c);
	//Vector / Scalar
	const Point3f operator/(const float& c) const;
	Point3f& operator/=(const float& c);
	//Vector * Vector (dot product)
	const Point3f operator*(const Point3f& v) const;
	Point3f& operator*=(const Point3f& v);
	//Vector ^ Vector (cross product)
	const Point3f operator^(const Point3f& v) const;
	Point3f& operator^=(const Point3f& v);
	//Add and subtract scalar
	const Point3f operator+(const float& c) const;
	Point3f& operator+=(const float& c);
	const Point3f operator-(const float& c) const;
	Point3f& operator-=(const float& c);
	//Add and subtract vector
	const Point3f operator+(const Point3f& v) const;
	Point3f& operator+=(const Point3f& v);
	const Point3f operator-(const Point3f& v) const;
	Point3f& operator-=(const Point3f& v);
	//Unary negation
	const Point3f operator-() const;
	//Array subscripts
	float& operator[](int idx);
	const float& operator[](int idx) const;
	//Equality
	bool operator==(const Point3f& v) const;
	bool operator!=(const Point3f& v) const;
	bool operator<(const Point3f& v) const;
	bool operator>(const Point3f& v) const;
	bool operator<=(const Point3f& v) const;
	bool operator>=(const Point3f& v) const;
	
	//Operations
	void setPosition(float val);
	void setPosition(float x, float y, float z);
	void setPosition(const Point3f &v);
};

//Scalar * vector
const Point3f operator*(const float& c, const Point3f& v);
const Point3f operator/(const float& c, const Point3f& v);

#endif

