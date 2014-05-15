#ifndef MATRIX2F_H
#define	MATRIX2F_H

class Matrix2f {
public:
	//[column][row] better for cpu caching?
	float data[2][2];
	
	Matrix2f(float i11, float i12, float i21, float i22);
	Matrix2f(const Matrix2f& orig);
	virtual ~Matrix2f();
	
	float determinant() const;
	const Matrix2f transpose() const;
	const Matrix2f inverse() const;
	//Singular value decomposition, where this = w.diag_product(e)*v.transpose()
	void svd(Matrix2f* w, Vector2f* e, Matrix2f* v) const;
	//Matrix * Matrix, where v represents a diagonal matrix
	void diag_product(const Vector2f& v);
	
	//Matrix / Scalar
	const Matrix2f operator/(const float& c) const;
	Matrix2f& operator/=(const float& c);
	//Matrix * Matrix
	const Matrix2f operator*(const Matrix2f& m) const;
	//Matrix * Vector
	const Vector2f operator*(const Vector2f& v) const;
};

const Matrix2f operator/(const float& c, const Matrix2f& m);

#endif

