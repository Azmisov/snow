#ifndef MAIN_H
#define	MAIN_H

#include "glfw3/glfw3.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include "Particle.h"
#include "Collider.h"
#include "PointCloud.h"
#include "Grid.h"

#define WIN_SIZE 640
#define WIN_METERS 1.5

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void redraw();
void *simulate(void *args);

#endif

