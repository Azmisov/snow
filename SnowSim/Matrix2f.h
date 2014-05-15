#ifndef MATRIX2F_H
#define	MATRIX2F_H

#include <cstring>
#include "Vector2f.h"

class Matrix2f {
public:
	//[column][row] better for cpu caching?
	float data[2][2];
	
	Matrix2f();
	Matrix2f(float i11, float i12, float i21, float i22);
	Matrix2f(const Matrix2f& m);
	Matrix2f(float data[2][2]);
	virtual ~Matrix2f();
	static Matrix2f identity(){
		return Matrix2f(1, 0, 0, 1);
	}
	
	void setData(const Matrix2f& m);
	void setData(float data[2][2]);
	void setData(const float data[2][2]);
	void setData(float i11, float i12, float i21, float i22);
	
	float determinant() const;
	const Matrix2f transpose() const;
	const Matrix2f inverse() const;
	//Singular value decomposition, where this = w.diag_product(e)*v.transpose()
	void svd(Matrix2f* w, Vector2f* e, Matrix2f* v) const;
	
	//DIAGONAL MATRIX OPERATIONS
	//Matrix * Matrix
	void diag_product(const Vector2f& v);
	//Matrix * Matrix^-1
	void diag_product_inv(const Vector2f& v);
	//Matrix - Matrix
	void diag_difference(const float& c);
	void diag_difference(const Vector2f& v);
	//Matrix + Matrix
	void diag_sum(const float& c);
	void diag_sum(const Vector2f& v);
	
	//SCALAR OVERLOADS
	//Matrix / Scalar
	const Matrix2f operator/(const float& c) const;
	Matrix2f& operator/=(const float& c);
	//Matrix * Scalar
	const Matrix2f operator*(const float& c) const;
	Matrix2f& operator*=(const float& c);
	//Matrix - Scalar
	const Matrix2f operator-(const float& c) const;
	Matrix2f& operator-=(const float& c);
	//Matrix + Scalar
	const Matrix2f operator+(const float& c) const;
	Matrix2f& operator+=(const float& c);
	
	//VECTOR OVERLOADS
	//Matrix + Matrix
	const Matrix2f operator+(const Matrix2f& m) const;
	Matrix2f& operator+=(const Matrix2f& m);
	//Matrix - Matrix
	const Matrix2f operator-(const Matrix2f& m) const;
	Matrix2f& operator-=(const Matrix2f& m);
	//Matrix * Matrix
	const Matrix2f operator*(const Matrix2f& m) const;
	//Matrix * Vector
	const Vector2f operator*(const Vector2f& v) const;
};

const Matrix2f operator/(const float& c, const Matrix2f& m);
const Matrix2f operator*(const float& c, const Matrix2f& m);
const Matrix2f operator+(const float& c, const Matrix2f& m);
const Matrix2f operator-(const float& c, const Matrix2f& m);

#endif

