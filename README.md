# PANG_model
## Aeroelastic modelling tool with 'Polynomial Approximated Nonlinear Geometry' (PANG)

Nonlinear aeroelastic modelling tool developed for S. Jayatilake, M. Lowenberg, B. Woods, B. Titurus, 'Nonlinear aeroelastic modelling and analysis of a geometrically nonlinear wing with combined unsteady sectional and lifting line aerodynamics', Nonlinear Dynamics, 2025. (https://link.springer.com/article/10.1007/s11071-025-10936-4). 

The current version includes the essential functionalities to build a model and call the aeroelastic systems. Two example cases are provided:

**example1**: modelling a 777 like wing with user-defined problem parameters,

**example2_BAFF**: shows an example with a BAFF aircraft file as an input for the problem parameters.

In each case, three example scripts are included:

`example_build.m`: this shows the process for building a model with some example parameters and generating and saving analysis modules.

`demo_ONERA.m`: this is an example of using the embedded aerodynamic model to perform some analyses.

`demoRun_extForce.m` this is an example for using the module for an external aerodynamic models, including the functionality for calling the deformed aerodyanamic meshes.

## Getting started

Add the `tbx` folder to the MATLAB path.

## Model building

Follow the example script `example_build.m`. This process is centralised around the `buildBase.m` class in the `buildSystem` namespace folder (accessed using `obj = buildSystem.buildBase`). 

## Analysis modules

Two analysis modules are available for aeroelastic work.

`run1 = analysis.oneraBase(obj)` allows you to generate an instance of this class (model-specific) to implement analysis using the embedded unsteady aerodynamic model, with finite span effects.

`run2 = analysis.extAeroBase(obj)` ...same for implementing an external aerodynamic model. This require the user to provide spanwise follower lift, drag forces and the moments. There're functionalities for getting the deformed aerodynamic meshes for external computations.
