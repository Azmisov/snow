#include "glfw3/glfw3.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define WIN_W 640
#define WIN_H 480

using namespace std;

static void error_callback(int, const char*);
void key_callback(GLFWwindow*, int, int, int, int);
void redraw();
void *simulate(void *args);

//Old and new time values for each timestep
double old_time, new_time = glfwGetTime();
bool dirty_buffer = true;

/*
 * 
 */
int main(int argc, char** argv) {
	//Create GLFW window
	GLFWwindow* window;
	glfwSetErrorCallback(error_callback);
	if (!glfwInit())
		exit(EXIT_FAILURE);
	window = glfwCreateWindow(WIN_W, WIN_H, "Snow Simulator", NULL, NULL);
	if (!window){
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	
	//Center window on screen
	const GLFWvidmode* monitor = glfwGetVideoMode(glfwGetPrimaryMonitor());
	glfwSetWindowPos(window, (monitor->width-WIN_W)/2, (monitor->height-WIN_H)/2);
	
	//Setup OpenGL context
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, WIN_W, WIN_H);
	
	//Create default simulation loop
	pthread_t sim_thread;
	pthread_create(&sim_thread, NULL, simulate, NULL);
	
	//Drawing & event loop
	while (!glfwWindowShouldClose(window)){
		if (dirty_buffer)
			redraw();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	//Exit
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
	return 0;
}
//Print all errors to console
static void error_callback(int error, const char* description){
	printf("\nError: %s",description);
}
//Key listener
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
    switch (key){
		default: return;
	}
}

/**
 * Draw the scene
 */
void redraw(){
	
}

/**
 * Computations at each time step
 */
void *simulate(void *args){
	cout << "Starting simulation thread..." << endl;
	
	//simulation variables
	int DIMS = 3;					//spacial dimensions to solve in
	double bounds[DIMS*2];			//bounding volume [xmin, xmax, ymin, ymax, ...]
	
	//initialize the domain
	
	
	pthread_exit(NULL);
}