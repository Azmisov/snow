#ifndef MAIN_H
#define	MAIN_H

#include "glfw3/glfw3.h"
#include "freeimage/FreeImage.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>
#include "Particle.h"
#include "PointCloud.h"
#include "Grid.h"
#include "SimConstants.h"
#include "Shape.h"

float TIMESTEP;

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void mouse_callback(GLFWwindow*, int, int, int);
void redraw();
void start_simulation();
void *simulate(void *args);
float adaptive_timestep();
void save_buffer(int time);

//Shape stuff
void create_new_shape();
void remove_all_shapes();

#endif

