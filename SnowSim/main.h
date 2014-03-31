#ifndef MAIN_H
#define	MAIN_H

#include "glfw3/glfw3.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <vector>
#include "Particle.h"
#include "Collider.h"

#define WIN_W 640
#define WIN_H 480

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void redraw();
/**
 * Create a cube point cloud
 * @param size number of particles along each axis
 * @param spacing spacing between particles
 * @return point cloud
 */
std::vector<Particle>* createCube(int size, float spacing);
/**
 * Computations at each time step
 */
void *simulate(void *args);

#endif

