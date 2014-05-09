#include "Point3f.h"
#include <cmath>
#include <stdio.h>
#include "float.h"

//CONSTRUCTORS
Point3f::Point3f(){
	this->setPosition(0, 0, 0);
}
Point3f::Point3f(float val){
	this->setPosition(val);
}
Point3f::Point3f(float x, float y, float z){
	this->setPosition(x, y, z);
}
Point3f::Point3f(const Point3f& orig) {
	this->setPosition(orig[0], orig[1], orig[2]);
}
Point3f::~Point3f(){
	
}

void Point3f::normalize(){
	double sum = this->length();
	for (int i=0; i<3; i++)
		loc[i] /= sum;
}
float Point3f::dot(const Point3f &v) const{
	return v.loc[0]*loc[0] + v.loc[1]*loc[1] + v.loc[2]*loc[2];
}
float Point3f::sum() const{
	return loc[0] + loc[1] + loc[2];
}
float Point3f::length() const{
	double sum = 0;
	for (int i=0; i<3; i++)
		sum += loc[i]*loc[i];
	return sqrt(sum);
}

//OVERLOADS
const Point3f operator*(const float& c, const Point3f& v){
	return Point3f(v)*c;
}
const Point3f operator/(const float& c, const Point3f& v){
	return Point3f(v)/c;
}
//Vector * Scalar
const Point3f Point3f::operator*(const float& c) const{
	return Point3f(*this) *= c;
}
Point3f& Point3f::operator*=(const float& c){
	loc[0] *= c;
	loc[1] *= c;
	loc[2] *= c;
	return *this;
}
//Vector / Scalar
const Point3f Point3f::operator/(const float& c) const{
	return Point3f(*this) /= c;
}
Point3f& Point3f::operator/=(const float& c){
	loc[0] /= c;
	loc[1] /= c;
	loc[2] /= c;
	return *this;
}
//Vector * Vector (dot product)
const Point3f Point3f::operator*(const Point3f& v) const{
	return Point3f(*this) *= v;
}
Point3f& Point3f::operator*=(const Point3f& v){
	loc[0] *= v.loc[0];
	loc[1] *= v.loc[1];
	loc[2] *= v.loc[2];
	return *this;
}
//Vector ^ Vector (cross product)
const Point3f Point3f::operator^(const Point3f& v) const{
	return Point3f(*this) ^= v;
}
Point3f& Point3f::operator^=(const Point3f& v){
	float v1 = loc[1]*v.loc[2] - loc[2]*v.loc[1],
		  v2 = loc[2]*v.loc[0] - loc[0]*v.loc[2],
		  v3 = loc[0]*v.loc[1] - loc[1]*v.loc[0];
	loc[0] = v1;
	loc[1] = v2;
	loc[2] = v3;
	return *this;
}
//Add and subtract
const Point3f Point3f::operator+(const float& c) const{
	return Point3f(*this) += c;
}
Point3f& Point3f::operator+=(const float& c){
	loc[0] += c;
	loc[1] += c;
	loc[2] += c;
	return *this;
}
const Point3f Point3f::operator-(const float& c) const{
	return Point3f(*this) -= c;
}
Point3f& Point3f::operator-=(const float& c){
	loc[0] -= c;
	loc[1] -= c;
	loc[2] -= c;
	return *this;
}
const Point3f Point3f::operator+(const Point3f& v) const{
	return Point3f(*this) += v;
}
Point3f& Point3f::operator+=(const Point3f& v){
	loc[0] += v.loc[0];
	loc[1] += v.loc[1];
	loc[2] += v.loc[2];
	return *this;
}
const Point3f Point3f::operator-(const Point3f& v) const{
	return Point3f(*this) -= v;
}
Point3f& Point3f::operator-=(const Point3f& v){
	loc[0] -= v.loc[0];
	loc[1] -= v.loc[1];
	loc[2] -= v.loc[2];
	return *this;
}
//Unary negation
const Point3f Point3f::operator-() const{
	return Point3f(-loc[0], -loc[1], -loc[2]);
}
//Array subscripts
float& Point3f::operator[](int idx){
	return loc[idx];
}
const float& Point3f::operator[](int idx) const{
	return loc[idx];
}
//Comparisons, since floating point is a mess
bool Point3f::operator==(const Point3f& v) const{
	return fabs(loc[0] - v.loc[0]) < FLT_EPSILON &&
		   fabs(loc[1] - v.loc[1]) < FLT_EPSILON &&
		   fabs(loc[2] - v.loc[2]) < FLT_EPSILON;
}
bool Point3f::operator!=(const Point3f& v) const{
	return !(*this == v);
}
bool Point3f::operator<(const Point3f& v) const{
	return loc[0] + FLT_EPSILON < v.loc[0] &&
		   loc[1] + FLT_EPSILON < v.loc[1] &&
		   loc[2] + FLT_EPSILON < v.loc[2];
}
bool Point3f::operator>(const Point3f& v) const{
	return v < *this;
}
bool Point3f::operator<=(const Point3f& v) const{
	return !(*this > v);
}
bool Point3f::operator>=(const Point3f& v) const{
	return !(*this < v);
}

//OPERATIONS
void Point3f::setPosition(float val){
	setPosition(val, val, val);
}
void Point3f::setPosition(float x, float y, float z){
	loc[0] = x;
	loc[1] = y;
	loc[2] = z;
}
void Point3f::setPosition(const Point3f &v){
	setPosition(v.loc[0], v.loc[1], v.loc[2]);
}