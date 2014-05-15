#include "Matrix2f.h"

Matrix2f::Matrix2f(){
	memset(data, 0, sizeof(float)*4);
}
Matrix2f::Matrix2f(float i11, float i12, float i21, float i22){
	setData(i11, i12, i21, i22);
}
Matrix2f::Matrix2f(const Matrix2f& orig){
	setData(orig);
}
Matrix2f::Matrix2f(float data[2][2]){
	setData(data);
}
Matrix2f::~Matrix2f(){}

void Matrix2f::setData(const Matrix2f& m){
	setData(m.data);
}
void Matrix2f::setData(float data[2][2]){
	memcpy(data, data, sizeof(float)*4);
}
void Matrix2f::setData(const float data[2][2]){
	setData(data[0][0], data[1][0], data[0][1], data[1][1]);
}
void Matrix2f::setData(float val){
	setData(val, val, val, val);
}
void Matrix2f::setData(float i11, float i12, float i21, float i22){
	data[0][0] = i11;
	data[0][1] = i21;
	data[1][0] = i12;
	data[1][1] = i22;
}

const float Matrix2f::determinant() const{
	return data[0][0]*data[1][1] - data[0][1]*data[1][0];
}
const Matrix2f Matrix2f::transpose() const{
	return Matrix2f(data[0][0], data[0][1], data[1][0], data[1][1]);
}
const Matrix2f Matrix2f::inverse() const{
	float det = determinant();
	return Matrix2f(
		data[1][1]/det,
		-data[1][0]/det,
		-data[0][1]/det,
		data[0][0]/det
	);
}
void Matrix2f::svd(Matrix2f* w, Vector2f* e, Matrix2f* v) const{
	//TODO: compute singular value decomposition
}

//DIAGONAL MATRIX OPERATIONS
//Matrix * Matrix
void Matrix2f::diag_product(const Vector2f& v){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] *= v[i];
	}
}
//Matrix * Matrix^-1
void Matrix2f::diag_product_inv(const Vector2f& v){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] /= v[i];
	}
}
//Matrix - Matrix
void Matrix2f::diag_difference(const float& c){
	for (int i=0; i<2; i++)
		data[i][i] -= c;
}
void Matrix2f::diag_difference(const Vector2f& v){
	for (int i=0; i<2; i++)
		data[i][i] -= v[i];
}
//Matrix + Matrix
void Matrix2f::diag_sum(const float& c){
	for (int i=0; i<2; i++)
		data[i][i] += c;
}
void Matrix2f::diag_sum(const Vector2f& v){
	for (int i=0; i<2; i++)
		data[i][i] += v[i];
}

//SCALAR OVERLOADS
//Matrix / Scalar
const Matrix2f operator/(const float& c, const Matrix2f& m){
	return Matrix2f(m)/c;
}
const Matrix2f Matrix2f::operator/(const float& c) const{
	return Matrix2f(*this) /= c;
}
Matrix2f& Matrix2f::operator/=(const float& c){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] /= c;
	}
	return *this;
}
//Matrix * Scalar
const Matrix2f operator*(const float& c, const Matrix2f& m){
	return Matrix2f(m)*c;
}
const Matrix2f Matrix2f::operator*(const float& c) const{
	return Matrix2f(*this) *= c;
}
Matrix2f& Matrix2f::operator*=(const float& c){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] *= c;
	}
	return *this;
}
//Matrix + Scalar
const Matrix2f operator+(const float& c, const Matrix2f& m){
	return Matrix2f(m)+c;
}
const Matrix2f Matrix2f::operator+(const float& c) const{
	return Matrix2f(*this) += c;
}
Matrix2f& Matrix2f::operator+=(const float& c){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] += c;
	}
	return *this;
}
//Matrix - Scalar
const Matrix2f operator-(const float& c, const Matrix2f& m){
	return Matrix2f(m)-c;
}
const Matrix2f Matrix2f::operator-(const float& c) const{
	return Matrix2f(*this) -= c;
}
Matrix2f& Matrix2f::operator-=(const float& c){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] -= c;
	}
	return *this;
}

//VECTOR OVERLOADS
//Matrix + Matrix
const Matrix2f Matrix2f::operator+(const Matrix2f& m) const{
	return Matrix2f(*this) += m;
}
Matrix2f& Matrix2f::operator+=(const Matrix2f& m){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] += m.data[i][j];
	}
	return *this;
}
//Matrix - Matrix
const Matrix2f Matrix2f::operator-(const Matrix2f& m) const{
	return Matrix2f(*this) -= m;
}
Matrix2f& Matrix2f::operator-=(const Matrix2f& m){
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++)
			data[i][j] -= m.data[i][j];
	}
	return *this;
}
//Matrix * Matrix
const Matrix2f Matrix2f::operator*(const Matrix2f& m) const{
	Matrix2f out = Matrix2f(*this);
	//Columns of m.data
	for (int i=0; i<2; i++){
		//Rows of data
		for (int j=0; j<2; j++){
			//Individual entries of each
			out.data[i][j] = data[0][j]*m.data[i][0];
			for (int k=1; k<2; k++)
				out.data[i][j] += data[k][j]*m.data[i][k];
		}
	}
	return out;
}
//Matrix * Vector
const Vector2f Matrix2f::operator*(const Vector2f& v) const{
	return Vector2f(
		data[0][0]*v[0] + data[1][0]*v[1],
		data[0][1]*v[0] + data[1][1]*v[1]
	);
}
