# DEMPgen
Event generator for Deep Exclusive Meson Production



## Building

To build the event generator, the compiler needs acces to both a compiled installation of CERN ROOT and its source code. ROOT version 6.08.06 or later is supported, and must be installed with the MathMore package enabled. Be sure not to change the location of either the ROOT source or compiled files after installation, as this will interfere with ROOT's built in CMake configurators.

CMake is also required. CMake 2.8 is the minimum supported version, and CMake 3 has been tested as well.

After downloading the source create a build directory and cd to it. Take note of the location of the source directory (where CMakeLists.txt should be stored) and run the commands:

mdkir build.
cd build
cmake ..
make

The event generator can now be run using the following command in the data/ directory.

cd data/
../build/main ../Config.json


## Configuration

The file Config.json contains all the configuration options. Use this as a template for other configuration files, which may be given as an argument to the event generator.

## Output

The event data is output to the configured file location relative to the build directory. The TTree in this file contains all kinematic data for all particles in the laboratory rest frame. Variables with the prefix "Vert" represent values read at the interaction vertex. Values with the prefix "Lab" represent values read after all correcting effects (multiple scattering, ionization, etc.) have been applied. 

## Sources and Classes

### Particle Class

The Particle class inherits from the root TLorentzVector class so that its functions are immediately accessible through the Particle objects. It's member variables include its mass, charge, and its GEMC compatible id.

### DEMPEvent Class

The Event class contains all particles of a single event, in one frame of reference. It includes methods enabling transformation to other frames of reference and coordinate systems.

### Asymmetry Class

The Asymmetry class calculates asymmetry amplitudes based on Monte-Carlo data by Goloskokov and Kroll. The asymmetry objects are persistent between events.

### Customrand Class

A class to hold various distribution functions for randomly generated variables, such as scattered electron energy. These functions should be persistent between events, and will provide fast random numbers.

### ScatteredParticleGen

Stores the kinematic ranges for the scattered electron and generates them with random energy and direction within the range, using sphere point picking.

### TargetGen

Generates the target neutron with Fermi momentum (when enabled) and proton for FSI (when enabled, may also have Fermi momentum).

### ProductGen

This class reads in the kinematic variables for the incident, target, and scattered particles and uses conservation laws to solve for the remaining two particles. The pion is first given a random direction. A fast root finding algorithm then calculates the pion's momentum magnitude, which is then used to find the proton's momentum. Pion direction may also be passed as an argument for debugging purposes.

### SigmaCalc

This class returns the cross sections and weights for the event. It acts as an interface between the current version of the event generator and header files from the old event generator (seen under branch "original"), which contain the parameterization of the cross sections.

### TreeBuilder

Manages the ROOT TTree object to be stored in the output .root file. Contains methods to easily add all kinematics data stored in a particle or DEMPEvent object. 

### Matter Effects

Transforms particle kinematics based on three effects caused by transition through matter: Ionization, Bremsstrahlung, and Multiple Scattering. 

### FSI

Computes the effect of final state interaction between the produced pion and one of the two protons of the target nucleus. The momentum of the outgoing pion and recoiled proton, as well as  the cross section of the interaction is calculated based on elastic pion-nucleon scattering.

## Acknowledgments

### JsonCpp

This project uses [JsonCpp](https://github.com/open-source-parsers/jsoncpp "JsonCpp Github") to read in configuration options. The amalgamated sources for JsonCpp are redistributed with this project in compliance with the MIT license.
