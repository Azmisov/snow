#ifndef MAIN_H
#define	MAIN_H

#include "glfw3/glfw3.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include "Particle.h"
#include "Collider.h"
#include "PointCloud.h"
#include "Grid.h"

//Various compiler options
#define SUPPORTS_POINT_SMOOTH 0
#define SCREENCAST 0
#if SCREENCAST
#include <FreeImage.h>
#endif

#define WIN_SIZE 640
#define WIN_METERS 1.5

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void redraw();
void *simulate(void *args);
void save_buffer(int time);

#endif

