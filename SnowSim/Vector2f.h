#ifndef VECTOR2F_H
#define	VECTOR2F_H

#include "Matrix2f.h"

class Matrix2f;

class Vector2f {
public:
	//Variables
	float data[2];
	
	//Constructors
	Vector2f();
	Vector2f(float val);
	Vector2f(float x, float y);
	Vector2f(const Vector2f& orig);
	virtual ~Vector2f();
	
	//Operations
	void setData(float val);
	void setData(float x, float y);
	void setData(const Vector2f &v);
	
	void normalize();
	const float dot(const Vector2f &v) const;
	const float sum() const;
	const float product() const;
	const float length() const;
	const float length_squared() const;
	//Vector * Vector^T
	const Matrix2f outer_product(const Vector2f& v) const;
	
	//OVERLOADS
	//Unary negation
	const Vector2f operator-() const;
	//Array subscripts
	float& operator[](int idx);
	const float& operator[](int idx) const;
	
	//SCALAR OVERLOADS
	//Vector * Scalar
	const Vector2f operator*(const float& c) const;
	Vector2f& operator*=(const float& c);
	//Vector / Scalar
	const Vector2f operator/(const float& c) const;
	Vector2f& operator/=(const float& c);
	//Vector + Scalar
	const Vector2f operator+(const float& c) const;
	Vector2f& operator+=(const float& c);
	//Vector - Scalar
	const Vector2f operator-(const float& c) const;
	Vector2f& operator-=(const float& c);
	
	//VECTOR OVERLOADS
	//Vector / Vector (piecewise division)
	const Vector2f operator/(const Vector2f& v) const;
	Vector2f& operator/=(const Vector2f& v);
	//Vector * Vector (dot product)
	const Vector2f operator*(const Vector2f& v) const;
	Vector2f& operator*=(const Vector2f& v);
	//Vector ^ Vector (cross product)
	const Vector2f operator^(const Vector2f& v) const;
	Vector2f& operator^=(const Vector2f& v);
	//Vector + Vector
	const Vector2f operator+(const Vector2f& v) const;
	Vector2f& operator+=(const Vector2f& v);
	//Vector - Vector
	const Vector2f operator-(const Vector2f& v) const;
	Vector2f& operator-=(const Vector2f& v);
};

//Scalar operations
const Vector2f operator*(const float& c, const Vector2f& v);
const Vector2f operator/(const float& c, const Vector2f& v);
const Vector2f operator-(const float& c, const Vector2f& v);
const Vector2f operator+(const float& c, const Vector2f& v);

#endif

