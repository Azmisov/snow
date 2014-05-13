/**
 * A 2-dimensional data value
 */

#ifndef VECTOR2F_H
#define	VECTOR2F_H

class Vector2f {
public:
	//Variables
	float loc[2];
	
	//Constructors
	Vector2f();
	Vector2f(float val);
	Vector2f(float x, float y);
	Vector2f(const Vector2f& orig);
	virtual ~Vector2f();
	
	void normalize();
	float dot(const Vector2f &v) const;
	float sum() const;
	float product() const;
	float length() const;
	
	//Operator Overloads
	//Vector * Scalar
	const Vector2f operator*(const float& c) const;
	Vector2f& operator*=(const float& c);
	//Vector / Scalar
	const Vector2f operator/(const float& c) const;
	Vector2f& operator/=(const float& c);
	//Vector / Vector (piecewise division)
	const Vector2f operator/(const Vector2f& v) const;
	Vector2f& operator/=(const Vector2f& v);
	//Vector * Vector (dot product)
	const Vector2f operator*(const Vector2f& v) const;
	Vector2f& operator*=(const Vector2f& v);
	//Vector ^ Vector (cross product)
	const Vector2f operator^(const Vector2f& v) const;
	Vector2f& operator^=(const Vector2f& v);
	//Add and subtract scalar
	const Vector2f operator+(const float& c) const;
	Vector2f& operator+=(const float& c);
	const Vector2f operator-(const float& c) const;
	Vector2f& operator-=(const float& c);
	//Add and subtract vector
	const Vector2f operator+(const Vector2f& v) const;
	Vector2f& operator+=(const Vector2f& v);
	const Vector2f operator-(const Vector2f& v) const;
	Vector2f& operator-=(const Vector2f& v);
	//Unary negation
	const Vector2f operator-() const;
	//Array subscripts
	float& operator[](int idx);
	const float& operator[](int idx) const;
	//Equality
	bool operator==(const Vector2f& v) const;
	bool operator!=(const Vector2f& v) const;
	bool operator<(const Vector2f& v) const;
	bool operator>(const Vector2f& v) const;
	bool operator<=(const Vector2f& v) const;
	bool operator>=(const Vector2f& v) const;
	
	//Operations
	void setPosition(float val);
	void setPosition(float x, float y);
	void setPosition(const Vector2f &v);
};

//Scalar * vector
const Vector2f operator*(const float& c, const Vector2f& v);
const Vector2f operator/(const float& c, const Vector2f& v);

#endif

