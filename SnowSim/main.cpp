#include "main.h"

using namespace std;

//Old and new time values for each timestep
double old_time, new_time = glfwGetTime();
bool dirty_buffer = true;
int frame_count = 0,
	bsize = 3*WIN_SIZE*WIN_SIZE;
unsigned char* img_buffer;

//Simulation data
int point_size;
PointCloud* snow;
Grid* grid;

int main(int argc, char** argv) {
	//Conjugate test
	Matrix2f A = Matrix2f(5, 6, 7, 8);
	Vector2f b = Vector2f(84, 116);
	Vector2f x = Vector2f(b);
	
	Vector2f r = b-A*x;
	Vector2f p = r;
	Vector2f Ap = A*p;
	Vector2f Ar = A*r;
	double rAr = r.dot(Ar);
	
	int iterations = 0;
	Vector2f error = Vector2f(100);
	while (error.length() > 0.0001){
		float alpha = rAr / (Ap.dot(Ap));
		error = alpha*p;
		x += error;
		//Update residual
		r -= alpha*Ap;
		Ar = A*r;
		//Compute new gradient
		double temp = r.dot(Ar);
		double beta = temp / rAr;
		rAr = temp;
		p = r + beta*p;
		Ap = Ar + beta*Ap;
		
		iterations++;
	}
	
	cout << "ITERS = " << iterations << endl;
	cout << "Error = " << error.length() << endl;
	cout << x[0] << ", " << x[1] << endl;
	Vector2f b_test = A*x;
	cout << "TARGET = " << b[0] << ", " << b[1] << endl;
	cout << "RESULT = " << b_test[0] << ", " << b_test[1] << endl;
	
	return 0;
	
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
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	//Setup simulation data
	const float mpdim = .25;		//meters in each dimension
	const int ppdim = 60;			//particle count for each dimension
	snow = PointCloud::createSquare(mpdim, ppdim);
	snow->translate(Vector2f(.75-mpdim/2*(ppdim != 1), 1));
	//Adjust visualization size to fill area
	point_size = WIN_SIZE/WIN_METERS*mpdim/ppdim-1;
	if (point_size < 1)
		point_size = 1;
	else if (point_size > 20)
		point_size = 20;
	if (!SUPPORTS_POINT_SMOOTH)
		point_size *= 3;
	
	grid = new Grid(Vector2f(0), Vector2f(WIN_METERS, WIN_METERS), Vector2f(100), snow);
	//We need to estimate particle volumes before we start
	grid->initializeMass();	
	grid->calculateVolumes();
	
	//Create default simulation loop
	pthread_t sim_thread;
	pthread_create(&sim_thread, NULL, simulate, NULL);
	
	//Drawing & event loop
	//Create directory to save buffers in
#if SCREENCAST
	mkdir("../screencast/",0777);
	FreeImage_Initialise();
	img_buffer = new unsigned char[bsize];
#endif
	
	while (!glfwWindowShouldClose(window)){
		if (dirty_buffer){
			redraw();
			dirty_buffer = false;
#if SCREENCAST
			save_buffer(frame_count++);
#endif
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	//Exit
#if SCREENCAST
	FreeImage_DeInitialise();
	delete[] img_buffer;
#endif	
	
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
	glPointSize(point_size);
	glBegin(GL_POINTS);
	for (int i=0; i<snow->size; i++){
		Particle& p = snow->particles[i];
		//We can use the particle's density to vary color
		//Max density set to 160+DENSITY
		float density = 1 - p.density/(160+DENSITY);
		glColor3f(0, density < 0 ? 0 : density, 1);
		glVertex2fv(p.position.loc);
	}
	glEnd();
	glDisable(GL_POINT_SMOOTH);
}

void *simulate(void *args){
	clock_t start = clock();
	cout << "Starting simulation..." << endl;
	Vector2f gravity = Vector2f(0, -9.8);

	struct timespec sleep_duration;
	sleep_duration.tv_sec = 0;
	sleep_duration.tv_nsec = TIMESTEP*1e9;
	
	while(true){
	//for (int i=0; i<550; i++){
		//Initialize FEM grid
		grid->initializeMass();
		grid->initializeVelocities();
		//Compute grid velocities
		grid->explicitVelocities(gravity);
		//Map back to particles
		grid->updateVelocities();
		//Update particle data
		snow->update();
		//Redraw snow
		dirty_buffer = true;
		//Delay... (if doing realtime visualization)
		//nanosleep(&sleep_duration, NULL);
	}

	cout << "Simulation complete: " << (clock()-start)/(float) CLOCKS_PER_SEC << " seconds" << endl;
	pthread_exit(NULL);
}

#if SCREENCAST
void save_buffer(int time){
	FILE *file;
	char fname[32];
	
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	sprintf(fname, "../screencast/t_%04d.png", time);
	printf("%s\n", fname);
	
	//Copy the image to buffer
	glReadBuffer(GL_BACK_LEFT);
	glReadPixels(0, 0, WIN_SIZE, WIN_SIZE, GL_BGR, GL_UNSIGNED_BYTE, img_buffer);
	FIBITMAP* img = FreeImage_ConvertFromRawBits(
		img_buffer, WIN_SIZE, WIN_SIZE, 3*WIN_SIZE,
		24, 0xFF0000, 0x00FF00, 0x0000FF, false
	);
	FreeImage_Save(FIF_PNG, img, fname, 0);
	FreeImage_Unload(img);
}
#endif
