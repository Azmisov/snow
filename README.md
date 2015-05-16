Snow Simulation
===============

This is an implementation of "A Material Point Method for Snow Simulation" (Stomakhin et al., 2013). There is a 2D version for OpenGL and a 3D version for Houdini.

[Demo video is on YouTube:](https://www.youtube.com/watch?v=13MqvmScNCc)

[![Demo video](https://img.youtube.com/vi/13MqvmScNCc/0.jpg)](https://www.youtube.com/watch?v=13MqvmScNCc)

## 2D Simulator

You may need to modify the project build settings for it to work.
Build dependencies:
- GLFW: http://www.glfw.org/
- FreeImage: http://freeimage.sourceforge.net/

Prebuilt versions of these libraries are included for debian based systems (x86-64).

### Controls
- **Click:** adds a point for a snow shape
- **Enter:** finishes previous snow shape and starts a new one
- **C:** allows you to create a circle shape (first click sets origin, second click sets radius)
- **F12:** converts snow shapes to particles and starts the simulation
- **ESC:** stops the simulation and removes all snow

## 3D Simulator

A Houdini digital asset, **ramshorn_fx_mpm_snow_otl_stable.otl** has been created for simulation and rendering setup. You'll need to install this otl as well as the snow solver node plugin.  Source code is in **SIM_SnowSolver.c**. Run **setup.sh** to build the plugin (Note: you may need to modify setup.sh to point to your houdini installation directory). See **tutorial.txt** and **tutorial.hipnc** for a basic setup. 
