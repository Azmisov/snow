#include "main.h"
#include "PointCloud.h"
#include "Grid.h"
#include <time.h>

using namespace std;

//Old and new time values for each timestep
double old_time, new_time = glfwGetTime();
bool dirty_buffer = true;

//Simulation data
int point_size;
PointCloud snow;
Grid* grid;

int main(int argc, char** argv) {	
	//Create GLFW window
	GLFWwindow* window;
	glfwSetErrorCallback(error_callback);
	if (!glfwInit())
		exit(EXIT_FAILURE);
	window = glfwCreateWindow(WIN_SIZE, WIN_SIZE, "Snow Simulator", NULL, NULL);
	if (!window){
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	
	//Center window on screen
	const GLFWvidmode* monitor = glfwGetVideoMode(glfwGetPrimaryMonitor());
	glfwSetWindowPos(window, (monitor->width-WIN_SIZE)/2, (monitor->height-WIN_SIZE)/2);
	
	//Setup OpenGL context
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, WIN_SIZE, WIN_SIZE);
	glOrtho(0, WIN_METERS, 0, WIN_METERS, 0, 1);
	
	//Set default visualization parameters
	glClearColor(1, 1, 1, 1);
	
	//Setup simulation data
	const float mpdim = .25;		//meters in each dimension
	const int ppdim = 1;		//particle count for each dimension
	snow = *PointCloud::createSquare(mpdim, ppdim);
	snow.translate(Vector2f(.75-mpdim/2*(ppdim != 1), 1));
	//Adjust visualization size to fill area
	point_size = WIN_SIZE/WIN_METERS*mpdim/ppdim-1;
	if (point_size < 1)
		point_size = 1;
	else if (point_size > 20)
		point_size = 20;
	
	grid = new Grid(Vector2f(0), Vector2f(WIN_METERS, WIN_METERS), Vector2f(60), &snow);
	grid->initialize();	
	grid->calculateVolumes();	//only for first iteration
	
	//Create default simulation loop
	pthread_t sim_thread;
	pthread_create(&sim_thread, NULL, simulate, NULL);
	
	//Drawing & event loop
	while (!glfwWindowShouldClose(window)){
		if (dirty_buffer){
			redraw();
			dirty_buffer = false;
		}
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

void redraw(){
	glClear(GL_COLOR_BUFFER_BIT);
		
	//Grid nodes
	glPointSize(1);
	glColor3f(1, 0, 0);
	glBegin(GL_POINTS);
	for (int i=0; i<grid->size[0]; i++){
		for (int j=0; j<grid->size[1]; j++)
			glVertex2fv((grid->origin+grid->cellsize*Vector2f(i, j)).loc);
	}
	glEnd();
	
	//Snow particles
	glEnable(GL_POINT_SMOOTH);
	glColor3f(0, 0, 1);
	glPointSize(point_size);
	glBegin(GL_POINTS);
	for (int i=0; i<snow.size; i++)
		glVertex2fv(snow.particles[i].position.loc);
	glEnd();
	glDisable(GL_POINT_SMOOTH);
}

void *simulate(void *args){
	clock_t start = clock();
	cout << "Starting simulation..." << endl;
	Vector2f gravity = Vector2f(0, -9.8);

	for (int i=0; i<1; i++){
		//Initialize FEM grid
		grid->initialize();
		//Compute grid velocities
		grid->calculateVelocities(gravity);
		//Map back to particles
		grid->updateVelocities();
		//Update particle data
		//snow.update();
		//Redraw snow
		dirty_buffer = true;
	}

	cout << "Simulation complete: " << (clock()-start)/CLOCKS_PER_SEC << " seconds" << endl;
	pthread_exit(NULL);
}
