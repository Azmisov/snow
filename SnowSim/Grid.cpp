#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* object){
	obj = object;
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	nodes_length = size.product();
	nodes = new GridNode[nodes_length];
	node_area = cellsize.product();
}
Grid::Grid(const Grid& orig){}
Grid::~Grid(){
	delete[] nodes;
}

//Maps mass and velocity to the grid
void Grid::initializeMass(){
	//Reset the grid
	//If the grid is sparsely filled, it may be better to reset individual nodes
	//Also, not all these variables need to be zeroed, so... yeah
	memset(nodes, 0, sizeof(GridNode)*nodes_length);
	
	//Map particle data to grid
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//Particle position to grid coordinates
		//This will give errors if the particle is outside the grid bounds
		p.grid_position = (p.position - origin)/cellsize;
		float ox = p.grid_position[0], oy = p.grid_position[1];
		
		
		//Shape function gives a blending radius of two;
		//so we do computations within a 2x2 square for each particle
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			//Y-dimension interpolation
			float y_pos = oy-y,
				wy = Grid::bspline(y_pos),
				dy = Grid::bsplineSlope(y_pos);
			
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				//X-dimension interpolation
				float x_pos = ox-x,
					wx = Grid::bspline(x_pos),
					dx = Grid::bsplineSlope(x_pos);
				
				//Final weight is dyadic product of weights in each dimension
				float weight = wx*wy;
				p.weights[idx] = weight;
				
				//Weight gradient is a vector of partial derivatives
				p.weight_gradient[idx].setData(dx*wy, wx*dy);				
				//I don't know why we need to do this... JT did it, doesn't appear in tech paper
				p.weight_gradient[idx] /= cellsize;
				
				//Interpolate mass
				nodes[(int) (y*size[0]+x)].mass += weight*p.mass;
			}
		}
	}
}
void Grid::initializeVelocities(){
	//We interpolate velocity after mass, to conserve momentum
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Interpolate velocity
					int n = (int) (y*size[0]+x);
					//We could also do a separate loop to divide by nodes[n].mass only once
					nodes[n].velocity += p.velocity * w * p.mass;
					nodes[n].active = true;
				}
			}
		}
	}
	for (int i=0; i<nodes_length; i++){
		GridNode &node = nodes[i];
		if (node.active)
			node.velocity /= node.mass;
	}
	collisionGrid();
}
//Maps volume from the grid to particles
//This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes() const{
	//Estimate each particles volume (for force calculations)
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		//First compute particle density
		p.density = 0;
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Node density is trivial
					p.density += w * nodes[(int) (y*size[0]+x)].mass;
				}
			}
		}
		p.density /= node_area;
		//Volume for each particle can be found from density
		p.volume = p.mass / p.density;
	}
}
//Calculate next timestep velocities for use in implicit integration
void Grid::explicitVelocities(const Vector2f& gravity){
	//First, compute the forces
	//We store force in velocity_new, since we're not using that variable at the moment
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//Solve for grid internal forces
		Matrix2f energy = p.energyDerivative();
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Weight the force onto nodes
					int n = (int) (y*size[0]+x);
					nodes[n].velocity_new += energy*p.weight_gradient[idx];
				}
			}
		}
	}
	
	//Now we have all grid forces, compute velocities (euler integration)
	for (int i=0; i<nodes_length; i++){
		GridNode &node = nodes[i];
		if (node.active)
			node.velocity_new = node.velocity + TIMESTEP*(gravity - node.velocity_new/node.mass);
	}
	collisionGrid();
}

#if ENABLE_IMPLICIT
//Solve linear system for implicit velocities
void Grid::implicitVelocities(){
	//With an explicit solution, we compute vf = vi + (f[n]/m)*dt
	//But for implicit, we use the force at the next timestep, f[n+1]
	//Stomakhin interpolates between the two, using IMPLICIT_RATIO
	//If we call v* the explicit vf, we can do some algebra and get
	//	v* = vf - IMPLICIT_RATION*dt*(df/m)
	//The problem is, df (change in force from n to n+1) depends on vf,
	//so we can't just compute it directly; instead, we use an iterative
	//method (conjugate residuals) to find what vf should be. We make an
	//initial guess of what vf should be (setting it to v*) and then
	//iteratively refine our guess until the error is small enough.
	
	//INITIALIZE LINEAR SOLVE
	for (int idx=0; idx<nodes_length; idx++){
		GridNode& n = nodes[idx];
		n.imp_active = n.active;
		if (n.imp_active){
			//recomputeImplicitForces will compute Er, given r
			//Initially, we want vf - E*vf; so we'll temporarily set r to vf
			n.r.setData(n.velocity_new);
			//Also set the error to 1
			n.err.setData(1);
		}
	}
	//As said before, we need to compute vf-E*vf as our initial "r" residual
	recomputeImplicitForces();
	for (int idx=0; idx<nodes_length; idx++){
		GridNode& n = nodes[idx];
		if (n.imp_active){
			n.r = n.velocity_new - n.Er;
			//p starts out equal to residual
			n.p = n.r;
			//cache r.dot(Er)
			n.rEr = n.r.dot(n.Er);
		}
	}
	//Since we updated r, we need to recompute Er
	recomputeImplicitForces();
	//Ep starts out the same as Er
	for (int idx=0; idx<nodes_length; idx++){
		GridNode& n = nodes[idx];
		if (n.imp_active)
			n.Ep = n.Er;
	}
	
	//LINEAR SOLVE
	for (int i=0; i<MAX_IMPLICIT_ITERS; i++){
		bool done = true;
		for (int idx=0; idx<nodes_length; idx++){
			GridNode& n = nodes[idx];
			//Only perform calculations on nodes that haven't been solved yet
			if (n.imp_active){
				//Alright, so we'll handle each node's solve separately
				//First thing to do is update our vf guess
				float div = n.Ep.dot(n.Ep);
				float alpha = n.rEr / div;
				n.err = alpha*n.p;
				//If the error is small enough, we're done
				float err = n.err.length();
				if (err < MAX_IMPLICIT_ERR || err > MIN_IMPLICIT_ERR || isnan(err)){
					n.imp_active = false;
					continue;
				}
				else done = false;
				//Update vf and residual
				n.velocity_new += n.err;
				n.r -= alpha*n.Ep;
			}
		}
		//If all the velocities converged, we're done
		if (done) break;
		//Otherwise we recompute Er, so we can compute our next guess
		recomputeImplicitForces();
		//Calculate the gradient for our next guess
		for (int idx=0; idx<nodes_length; idx++){
			GridNode& n = nodes[idx];
			if (n.imp_active){
				float temp = n.r.dot(n.Er);
				float beta = temp / n.rEr;
				n.rEr = temp;
				//Update p
				n.p *= beta;
				n.p += n.r;
				//Update Ep
				n.Ep *= beta;
				n.Ep += n.Er;
			}
		}
	}
}
void Grid::recomputeImplicitForces(){
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				GridNode& n = nodes[(int) (y*size[0]+x)];
				if (n.imp_active){
					//I don't think there is any way to cache intermediary
					//results for reuse with each iteration, unfortunately
					n.force += p.deltaForce(n.r, p.weight_gradient[idx]);
				}
			}
		}
	}
	
	//We have delta force for each node; to get Er, we use the following formula:
	//	r - IMPLICIT_RATIO*TIMESTEP*delta_force/mass
	for (int idx=0; idx<nodes_length; idx++){
		GridNode& n = nodes[idx];
		if (n.imp_active)
			n.Er = n.r - IMPLICIT_RATIO*TIMESTEP/n.mass*n.force;
	}
}
#endif

//Map grid velocities back to particles
void Grid::updateVelocities() const{
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//We calculate PIC and FLIP velocities separately
		Vector2f pic, flip = p.velocity;
		//Also keep track of velocity gradient
		Matrix2f& grad = p.velocity_gradient;
		grad.setData(0.0);
		//VISUALIZATION PURPOSES ONLY:
		//Recompute density
		p.density = 0;
		
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					GridNode &node = nodes[(int) (y*size[0]+x)];
					//Particle in cell
					pic += node.velocity_new*w;
					//Fluid implicit particle
					flip += (node.velocity_new - node.velocity)*w;
					//Velocity gradient
					grad += node.velocity_new.outer_product(p.weight_gradient[idx]);
					//VISUALIZATION ONLY: Update density
					p.density += w * node.mass;
				}
			}
		}
		//Final velocity is a linear combination of PIC and FLIP components
		p.velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
		//VISUALIZATION: Update density
		p.density /= node_area;
	}
	collisionParticles();
}

void Grid::collisionGrid(){
	Vector2f delta_scale = Vector2f(TIMESTEP);
	delta_scale /= cellsize;
	for (int y=0, idx=0; y<size[1]; y++){
		for (int x=0; x<size[0]; x++, idx++){
			//Get grid node (equivalent to (y*size[0] + x))
			GridNode &node = nodes[idx];
			//Check to see if this node needs to be computed
			if (node.active){
				//Collision response
				//TODO: make this work for arbitrary collision geometry
				Vector2f new_pos = node.velocity_new*delta_scale + Vector2f(x, y);
				//Left border, right border
				if (new_pos[0] < BSPLINE_RADIUS || new_pos[0] > size[0]-BSPLINE_RADIUS-1){
					node.velocity_new[0] = 0;
					node.velocity_new[1] *= STICKY;
				}
				//Bottom border, top border
				if (new_pos[1] < BSPLINE_RADIUS || new_pos[1] > size[1]-BSPLINE_RADIUS-1){
					node.velocity_new[0] *= STICKY;
					node.velocity_new[1] = 0;
				}
			}
		}
	}
}
void Grid::collisionParticles() const{
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		Vector2f new_pos = p.grid_position + TIMESTEP*p.velocity/cellsize;
		//Left border, right border
		if (new_pos[0] < BSPLINE_RADIUS-1 || new_pos[0] > size[0]-BSPLINE_RADIUS)
			p.velocity[0] = -STICKY*p.velocity[0];
		//Bottom border, top border
		if (new_pos[1] < BSPLINE_RADIUS-1 || new_pos[1] > size[1]-BSPLINE_RADIUS)
			p.velocity[1] = -STICKY*p.velocity[1];
	}
}