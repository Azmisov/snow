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
#include "PointCloud.h"
#include "Grid.h"

//Various compiler options
#define LIMIT_FPS 0
#define FPS 30
#define REDRAW_EVERY LIMIT_FPS ? 1/TIMESTEP/FPS : 1;
#define REALTIME_PLAYBACK 0
#define SUPPORTS_POINT_SMOOTH 0
#define SCREENCAST 0
#if SCREENCAST
#include <FreeImage.h>
#endif

#define WIN_SIZE 640
#define WIN_METERS 1.5

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void mouse_callback(GLFWwindow*, int, int, int);
void redraw();
void start_simulation();
void *simulate(void *args);
void save_buffer(int time);

//Shape stuff
void create_new_shape();
void remove_all_shapes();

#endif

