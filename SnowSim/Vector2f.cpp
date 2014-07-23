#include "Vector2f.h"
#include <cmath>
#include <stdio.h>

//CONSTRUCTORS
Vector2f::Vector2f(){
	this->setData(0, 0);
}
Vector2f::Vector2f(float val){
	this->setData(val);
}
Vector2f::Vector2f(float x, float y){
	this->setData(x, y);
}
Vector2f::Vector2f(const Vector2f& orig) {
	this->setData(orig[0], orig[1]);
}
Vector2f::~Vector2f(){}

//OPERATIONS
void Vector2f::setData(float val){
	setData(val, val);
}
void Vector2f::setData(float x, float y){
	data[0] = x;
	data[1] = y;
}
void Vector2f::setData(const Vector2f &v){
	setData(v.data[0], v.data[1]);
}

void Vector2f::normalize(){
	*this /= this->length();
}
const float Vector2f::dot(const Vector2f &v) const{
	return v.data[0]*data[0] + v.data[1]*data[1];
}
const float Vector2f::sum() const{
	return data[0] + data[1];
}
const float Vector2f::product() const{
	return data[0] * data[1];
}
const float Vector2f::length() const{
	return sqrt(length_squared());
}
const float Vector2f::length_squared() const{
	double sum = 0;
	for (int i=0; i<2; i++)
		sum += data[i]*data[i];
	return sum;
}
const Matrix2f Vector2f::outer_product(const Vector2f& v) const{
	return Matrix2f(
		data[0]*v.data[0], data[0]*v.data[1],
		data[1]*v.data[0], data[1]*v.data[1]
	);
}

//OVERLOADS
//Array subscripts
float& Vector2f::operator[](int idx){
	return data[idx];
}
const float& Vector2f::operator[](int idx) const{
	return data[idx];
}

//SCALAR OVERLOADS
//Vector * Scalar
const Vector2f operator*(const float& c, const Vector2f& v){
	return Vector2f(v)*c;
}
const Vector2f Vector2f::operator*(const float& c) const{
	return Vector2f(*this) *= c;
}
Vector2f& Vector2f::operator*=(const float& c){
	data[0] *= c;
	data[1] *= c;
	return *this;
}
//Vector / Scalar
const Vector2f operator/(const float& c, const Vector2f& v){
	return Vector2f(v)/c;
}
const Vector2f Vector2f::operator/(const float& c) const{
	return Vector2f(*this) /= c;
}
Vector2f& Vector2f::operator/=(const float& c){
	data[0] /= c;
	data[1] /= c;
	return *this;
}
//Vector + Scalar
const Vector2f operator+(const float& c, const Vector2f& v){
	return Vector2f(v)+c;
}
const Vector2f Vector2f::operator+(const float& c) const{
	return Vector2f(*this) += c;
}
Vector2f& Vector2f::operator+=(const float& c){
	data[0] += c;
	data[1] += c;
	return *this;
}
//Vector - Scalar
const Vector2f operator-(const float& c, const Vector2f& v){
	return Vector2f(v)-c;
}
const Vector2f Vector2f::operator-(const float& c) const{
	return Vector2f(*this) -= c;
}
Vector2f& Vector2f::operator-=(const float& c){
	data[0] -= c;
	data[1] -= c;
	return *this;
}

//VECTOR OVERLOADS
//Vector / Vector (piecewise division)
const Vector2f Vector2f::operator/(const Vector2f& v) const{
	return Vector2f(*this) /= v;
}
Vector2f& Vector2f::operator/=(const Vector2f& v){
	data[0] /= v.data[0];
	data[1] /= v.data[1];
	return *this;
}
//Vector * Vector (dot product)
const Vector2f Vector2f::operator*(const Vector2f& v) const{
	return Vector2f(*this) *= v;
}
Vector2f& Vector2f::operator*=(const Vector2f& v){
	data[0] *= v.data[0];
	data[1] *= v.data[1];
	return *this;
}
//Vector ^ Vector (cross product)
const Vector2f Vector2f::operator^(const Vector2f& v) const{
	return Vector2f(*this) ^= v;
}
Vector2f& Vector2f::operator^=(const Vector2f& v){
	//TODO: this may be incorrect...
	float v1 = data[0]*v.data[1],
		  v2 = -data[1]*v.data[0];
	data[0] = v1;
	data[1] = v2;
	return *this;
}
//Vector + Vector
const Vector2f Vector2f::operator+(const Vector2f& v) const{
	return Vector2f(*this) += v;
}
Vector2f& Vector2f::operator+=(const Vector2f& v){
	data[0] += v.data[0];
	data[1] += v.data[1];
	return *this;
}
//Vector - Vector
const Vector2f Vector2f::operator-(const Vector2f& v) const{
	return Vector2f(*this) -= v;
}
Vector2f& Vector2f::operator-=(const Vector2f& v){
	data[0] -= v.data[0];
	data[1] -= v.data[1];
	return *this;
}
//Unary negation
const Vector2f Vector2f::operator-() const{
	return Vector2f(-data[0], -data[1]);
}