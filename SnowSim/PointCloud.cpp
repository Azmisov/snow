#include "PointCloud.h"

PointCloud::PointCloud(){}
PointCloud::PointCloud(const PointCloud& orig){}
PointCloud::~PointCloud(){}

void PointCloud::scale(Vector2f origin, Vector2f scale){
	for (int i=0; i<size; i++){
		for (int j=0; j<2; j++){
			particles[i].position[j] = origin[j] + (particles[i].position[j]-origin[j])*scale[j];
		}
	}
}
void PointCloud::translate(Vector2f off){
	for (int i=0; i<size; i++){
		particles[i].position[0] += off[0];
		particles[i].position[1] += off[1];
	}
}
void PointCloud::bounds(float bounds[4]){
	bounds[0] = particles[0].position[0]; bounds[1] = bounds[0];
	bounds[2] = particles[0].position[1]; bounds[3] = bounds[2];
	for (int i=0; i<size; i++){
		Vector2f &p = particles[i].position;
		//X-bounds
		if (p[0] < bounds[0])
			bounds[0] = p[0];
		else if (p[0] > bounds[1])
			bounds[1] = p[0];
		//Y-bounds
		if (p[1] < bounds[2])
			bounds[2] = p[1];
		else if (p[1] > bounds[3])
			bounds[3] = p[1];
	}
}