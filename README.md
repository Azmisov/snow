Snow Simulation
===============

You may need to modify the project build settings for it to work.
Build dependencies:
- GLFW: http://www.glfw.org/
- FreeImage: http://freeimage.sourceforge.net/
Prebuilt versions of these libraries are included for Debian x86-64

This is an implementation of "A Material Point Method for Snow Simulation" (Stomakhin et al., 2013) in 2D.

Controls
========
- **Click:** adds a point for a snow shape
- **Enter:** finishes previous snow shape and starts a new one
- **C:** allows you to create a circle shape (first click sets origin, second click sets radius)
- **F12:** converts snow shapes to particles and starts the simulation
- **ESC:** stops the simulation and removes all snow
