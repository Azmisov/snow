#include "main.h"

using namespace std;

//Old and new time values for each timestep
double old_time, new_time = glfwGetTime();
bool dirty_buffer = true;
vector<Particle>* pcloud;

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
	
	//Setup simulation data
	pcloud = createCube(5, 1);
	
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

vector<Particle>* createCube(int size, float spacing){
	//point cloud definition; simple cube
	vector<Particle> *cloud = new vector<Particle>(size*size*size);
	for (int x=0, n=0; x<size; x++){
		for (int y=0; y<size; y++){
			for (int z=0; z<size; z++, n++){
				(*cloud)[n].setPosition(x*spacing, y*spacing, z*spacing);
			}
		}
	}
	return cloud;
}


void *simulate(void *args){
	cout << "Starting simulation thread..." << endl;
	
	//bounding volume [xmin, xmax, ymin, ymax, ...]
	double bounds[6];
	Particle* snow = &(*pcloud)[0];
	
	//initialize the domain
	int obj_size = pcloud->size();
	bounds[0] = snow[0].x; bounds[1] = snow[0].x;
	bounds[2] = snow[0].y; bounds[3] = snow[0].y;
	bounds[4] = snow[0].z; bounds[5] = snow[0].z;
	for (int i=0; i<obj_size; i++){
		Particle &p = snow[i];
		//X-bounds
		if (p.x < bounds[0])
			bounds[0] = p.x;
		else if (p.x > bounds[1])
			bounds[1] = p.x;
		//Y-bounds
		if (p.y < bounds[2])
			bounds[2] = p.y;
		else if (p.y > bounds[3])
			bounds[3] = p.y;
		//Z-bounds
		if (p.z < bounds[4])
			bounds[4] = p.z;
		else if (p.z > bounds[5])
			bounds[5] = p.z;
	}
	
	pthread_exit(NULL);
}