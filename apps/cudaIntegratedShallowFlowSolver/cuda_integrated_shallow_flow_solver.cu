// ======================================================================================
// Name                :    High-Performance Integrated Modelling System
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software. Legacy name: GeoClasses
// ======================================================================================
// Version             :    1.0.1
// Author              :    Xilin Xia
// Create Time         :    2014/10/04
// Update Time         :    2020/04/19
// ======================================================================================
// LICENCE: GPLv3 
// ======================================================================================


/*!
\file cuda_integrated_shallow_flow_solver.cu
\brief Source file for component test

*/


#include <iostream>
#include <cuda_runtime_api.h>
#include <cuda.h>

//These header files are the primitive types
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "cuda_arrays.h"
//These header files are for the fields
#include "mapped_field.h"
#include "cuda_mapped_field.h"
//These header files are for finite volume mesh
#include "mesh_fv_reduced.h"
#include "mesh_interface.h"
#include "cuda_mesh_fv.h"
//These header files are for input and output
#include "gisAsciiMesh_reader.h"
#include "gmsh_reader.h"
#include "field_reader.h"
#include "cuda_simple_writer.h"
#include "cuda_gisascii_writer.h"
//These header files are for shallow water equations advection
#include "cuda_advection_NSWEs.h"
//The header file for gradient
#include "cuda_gradient.h"
//The header file for hessian
#include "cuda_hessian.h"
//The header file for friction
#include "cuda_friction.h"
//The header file for limiter
#include "cuda_limiter.h"
//The header file for field algebra
#include "cuda_field_algebra.h"
//The header file for integrator
#include "cuda_integrators.h"
//The header file for time controllinh
#include "cuda_adaptive_time_control.h"
#include "cuda_device_query.h"

//using the name space for GeoClasses
using namespace GC;

int main(){

  deviceQuery();

  //creating an structured finite volume mesh, it is initialized from DEM
  std::shared_ptr<unstructuredFvMesh>  mesh = std::make_shared<unstructuredReducedFvMesh>(gisAsciiMeshReader("input/mesh/DEM.txt"));

  std::cout<<"Read in mesh successfully"<<std::endl;

  //creating mesh on device
  std::shared_ptr<cuUnstructuredFvMesh>  mesh_ptr_dev = std::make_shared<cuUnstructuredFvMesh>(fvMeshQueries(mesh));
  
  //Read in field data
  fvScalarFieldOnCell z_host(fvMeshQueries(mesh), completeFieldReader("input/field/", "z"));
  fvScalarFieldOnCell h_host(fvMeshQueries(mesh), completeFieldReader("input/field/", "h"));
  fvVectorFieldOnCell hU_host(fvMeshQueries(mesh), completeFieldReader("input/field/", "hU"));

  std::cout << "Read in field successfully" << std::endl;

  //h, z, hU
  cuFvMappedField<Scalar, on_cell> z(z_host,mesh_ptr_dev);
  cuFvMappedField<Scalar, on_cell> h(h_host,mesh_ptr_dev);
  cuFvMappedField<Vector, on_cell> hU(hU_host, mesh_ptr_dev);

  //Velocity
  cuFvMappedField<Vector, on_cell> u(hU, partial);

  //advections
  cuFvMappedField<Scalar, on_cell> h_advection(h, partial);
  cuFvMappedField<Vector, on_cell> hU_advection(hU, partial);

  //friction force
  cuFvMappedField<Vector, on_cell> friction_force(hU, partial);

  //flux
  cuFvMappedField<Scalar, on_halffacet> h_flux(mesh_ptr_dev);
  cuFvMappedField<Vector, on_halffacet> hU_flux(mesh_ptr_dev);

  //old values
  cuFvMappedField<Scalar, on_cell> h_old(h, partial);
  cuFvMappedField<Vector, on_cell> hU_old(hU, partial);
  cuFvMappedField<Scalar, on_cell> h_advection_old(h_advection, partial);
  cuFvMappedField<Vector, on_cell> hU_advection_old(hU_advection, partial);

  //gradients
  cuFvMappedField<Vector, on_cell> h_grad(mesh_ptr_dev);
  cuFvMappedField<Vector, on_cell> z_grad(mesh_ptr_dev);
  cuFvMappedField<Vector, on_cell> eta_grad(mesh_ptr_dev);
  cuFvMappedField<Tensor, on_cell> u_grad(mesh_ptr_dev);

  //modifying factor for gravity
  cuFvMappedField<Scalar, on_cell> phi_b(mesh_ptr_dev);


  //calculate topography gradient
  fv::cuGradient(z, z_grad);
  fv::cuGradientLimiterCartesian(z, z_grad);

  //hessian
  cuFvMappedField<Tensor, on_cell> z_hessian(mesh_ptr_dev);
  fv::cuHessianCartesian2D(z, z_hessian);

  //gravity
  cuFvMappedField<Scalar, on_cell> gravity(h, partial);
  //setting gravity to single value 9.81
  fv::cuUnaryOn(gravity, [] __device__ (Scalar& a) -> Scalar{return 9.81;}); 

  //centrifugal
  cuFvMappedField<Scalar, on_cell> centrifugal(h, partial);

  //gross vertical acceleration
  cuFvMappedField<Scalar, on_cell> gross_acc(h, partial);

  //surface elevation eta
  cuFvMappedField<Scalar, on_cell> eta(h, partial);

  Scalar dt_out = 0.5;
  Scalar t_out = 0.0;
  Scalar t_all = 0.0;
  Scalar CFL = 0.5;
  Flag rheology_type = 0;
  Flag gravity_correction_type = 0;
  bool curvature_on = true;
  bool filter_mass_flux = true;

  //creating uniform friction coefficients field
  cuFvMappedField<Scalar, on_cell> miu1(h, partial);
  cuFvMappedField<Scalar, on_cell> miu2(h, partial);
  cuFvMappedField<Scalar, on_cell> cohesion(h, partial);
  cuFvMappedField<Scalar, on_cell> rho(h, partial);

  Scalar _miu, _miu1, _miu2, _cohesion, _rho, _L, _beta, _U;

  //reading setup configuration
  std::ifstream setup;
  std::string line;
  std::string word;
  setup.open("input/setup.conf");
  if (!setup){
    std::cout << "error: unable to open input file: " << "input/setup.conf" << std::endl;
  }
  getline(setup, line);
  getline(setup, line);
  std::stringstream(line) >> t_all;
  getline(setup, line);
  getline(setup, line);
  std::stringstream(line) >> dt_out;
  getline(setup, line);
  getline(setup, line);
  std::stringstream(line) >> CFL;
  getline(setup, line);
  getline(setup, line);
  std::stringstream(line) >> rheology_type;
  switch (rheology_type){
  case 1:
    getline(setup, line);
    getline(setup, line);
    std::stringstream(line) >> _miu>> _cohesion>> _rho;
    if (_rho <= 1e-6){
      std::cout << "Fatal error: Density is 0!" << std::endl;
    }
    fv::cuUnaryOn(miu1, [=] __device__(Scalar& a) -> Scalar{ return _miu; });
    fv::cuUnaryOn(cohesion, [=] __device__(Scalar& a) -> Scalar{ return _cohesion; });
    fv::cuUnaryOn(rho, [=] __device__(Scalar& a) -> Scalar{ return _rho; });
    break;
  case 2:
    getline(setup, line);
    getline(setup, line);
    std::stringstream(line) >> _miu1 >> _miu2 >>_L >> _beta;
    fv::cuUnaryOn(miu1, [=] __device__(Scalar& a) -> Scalar{ return _miu1; });
    fv::cuUnaryOn(miu2, [=] __device__(Scalar& a) -> Scalar{ return _miu2; });
    break;
  case 3:
    getline(setup, line);
    getline(setup, line);
    std::stringstream(line) >> _miu1 >> _miu2 >> _U;
    fv::cuUnaryOn(miu1, [=] __device__(Scalar& a) -> Scalar{ return _miu1; });
    fv::cuUnaryOn(miu2, [=] __device__(Scalar& a) -> Scalar{ return _miu2; });
    break;
  }
  getline(setup, line);
  getline(setup, line);
  std::stringstream(line) >> gravity_correction_type;
  getline(setup, line);
  setup>>std::boolalpha >> curvature_on;
  getline(setup, line);
  getline(setup, line);
  setup >> std::boolalpha >> filter_mass_flux;
  getline(setup, line);
  setup.close();
  std::cout << "Read in configuration successfully" << std::endl;

  cuAdaptiveTimeControl2D time_controller(0.005, t_all, CFL);

  //ascii raster writer
  cuGisAsciiWriter raster_writer("input/mesh/DEM.txt");

  //write the initial profile
  cuSimpleWriterLowPrecision(z, "z", time_controller.current());
  cuSimpleWriterLowPrecision(h, "h", time_controller.current());
  cuSimpleWriterLowPrecision(hU, "hU", time_controller.current());
  cuSimpleWriterLowPrecision(eta, "elev", time_controller.current());

  double total_runtime = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  //print current time
  std::cout << time_controller.current() << std::endl;

  auto filter = [] __device__(Vector& a, Scalar& b) ->Vector{
    if (b <= 1e-6){
      return Vector(0.0);
    }
    else{
      return a;
    }
  };

  fv::cuBinaryOn(hU, h, filter);

  auto filter1 = [] __device__(Scalar& a) ->Scalar{
    if (a <= 1e-10){
      return 0.0;
    }
    else{
      return a;
    }
  };

  auto filter2 = [] __device__(Scalar& a, Vector& b) ->Scalar{
    if (sqrt(dot(b, b)) <= 1e-8){
      return 0.0;
    }
    else{
      return a;
    }
  };

  //Main loop
  do{

    cudaEventRecord(start);

    //First step of Runger-Kutta--------------------------------------

    //calculate the surface elevation
    fv::cuBinary(h, z, eta, [] __device__ (Scalar& a, Scalar& b) -> Scalar{return a + b;});

    //calculate the velocity
    auto divide = [] __device__ (Vector& a, Scalar& b) ->Vector{
      if(b >= 1e-10){
        return a/b;
      }else{
        return Vector(0.0);
      }
    };

    fv::cuBinary(hU, h, u, divide); 

    //store the old values, may need to be replaced by assining operators
    fv::cuBinaryOn(h_old, h, [] __device__ (Scalar& a, Scalar& b) -> Scalar{return b;});
    fv::cuBinaryOn(hU_old, hU, [] __device__ (Vector& a, Vector& b) -> Vector{return b;});

    //calculate gradient
    fv::cuGradient(h, h_grad);
    fv::cuGradient(eta, eta_grad);
    fv::cuGradient(u, u_grad);

    //Limiting gradient
    fv::cuGradientLimiterCartesian(h, h_grad);
    fv::cuGradientLimiterCartesian(eta, eta_grad);
    fv::cuGradientLimiterCartesian(u, u_grad);
    
    switch (gravity_correction_type){
    case 0:
      break;
    case 1:
      //setting gravity to single value 9.81
      fv::cuUnaryOn(gravity, [] __device__(Scalar& a) -> Scalar{ return 9.81; });
      fv::cuBinary(z_grad, gravity, gravity, [] __device__(Vector& a, Scalar& b) -> Scalar{ return b / (1.0 + dot(a, a)); });
      break;
    case 2:
      //calculating modifying factor
      fv::cuBinary(z_grad, eta_grad, phi_b, [] __device__(Vector& a, Vector& b) -> Scalar{ return fmax(0.0, fmin(1.0, 1 / (1.0 + dot(a, b)))); });
      //setting gravity to single value 9.81
      fv::cuUnaryOn(gravity, [] __device__(Scalar& a) -> Scalar{ return 9.81; });
      fv::cuBinary(phi_b, gravity, gravity, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return a*b; });
      break;
    }


    //calculate centrifugal force
    if (curvature_on){
      fv::cuBinary(u, z_hessian, centrifugal, [] __device__(Vector& a, Tensor& b) -> Scalar{ return dot(dot(b, a), a); });
      fv::cuBinary(z_grad, centrifugal, centrifugal, [] __device__(Vector& a, Scalar& b) -> Scalar{ return b / (1.0 + dot(a, a)); });
      fv::cuBinary(gravity, centrifugal, centrifugal, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return fmin(10.0*a, fmax(-0.9*a, b)); });
    }

    //calculate advection
    fv::cuAdvectionNSWEs2ndRobustCurv(gravity, centrifugal, h, z, u, h_grad, eta_grad, u_grad, h_advection, hU_advection);

    //multiply advection with -1
    fv::cuUnaryOn(h_advection, [] __device__ (Scalar& a) -> Scalar{return -1.0*a;});
    fv::cuUnaryOn(hU_advection, [] __device__ (Vector& a) -> Vector{return -1.0*a;});

    //integration
    fv::cuEulerIntegrator(h, h_advection, time_controller.dt(), time_controller.current());
    fv::cuEulerIntegrator(hU, hU_advection, time_controller.dt(), time_controller.current());

    //store the old values, may need to be replaced by assining operators
    fv::cuBinaryOn(h_advection_old, h_advection, [] __device__ (Scalar& a, Scalar& b) -> Scalar{return b;});
    fv::cuBinaryOn(hU_advection_old, hU_advection, [] __device__ (Vector& a, Vector& b) -> Vector{return b;});

    fv::cuUnaryOn(h, filter1);

    //Second step of Runger-Kutta--------------------------------------


    //calculate the surface elevation
    fv::cuBinary(h, z, eta, [] __device__ (Scalar& a, Scalar& b) -> Scalar{return a + b;});

    //calculate the velocity
    fv::cuBinary(hU, h, u, divide); 

    //calculate gradient
    fv::cuGradient(h, h_grad);
    fv::cuGradient(eta, eta_grad);
    fv::cuGradient(u, u_grad);

    //Limiting gradient
    fv::cuGradientLimiterCartesian(h, h_grad);
    fv::cuGradientLimiterCartesian(eta, eta_grad);
    fv::cuGradientLimiterCartesian(u, u_grad);

    switch (gravity_correction_type){
    case 0:
      break;
    case 1:
      //setting gravity to single value 9.81
      fv::cuUnaryOn(gravity, [] __device__(Scalar& a) -> Scalar{ return 9.81; });
      fv::cuBinary(z_grad, gravity, gravity, [] __device__(Vector& a, Scalar& b) -> Scalar{ return b / (1.0 + dot(a, a)); });
      break;
    case 2:
      //calculating modifying factor
      fv::cuBinary(z_grad, eta_grad, phi_b, [] __device__(Vector& a, Vector& b) -> Scalar{ return fmax(0.0, fmin(1.0, 1 / (1.0 + dot(a, b)))); });
      //setting gravity to single value 9.81
      fv::cuUnaryOn(gravity, [] __device__(Scalar& a) -> Scalar{ return 9.81; });
      fv::cuBinary(phi_b, gravity, gravity, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return a*b; });
      break;
    }

    if (curvature_on){
      //calculate centrifugal force
      fv::cuBinary(u, z_hessian, centrifugal, [] __device__(Vector& a, Tensor& b) -> Scalar{ return dot(dot(b, a), a); });
      fv::cuBinary(z_grad, centrifugal, centrifugal, [] __device__(Vector& a, Scalar& b) -> Scalar{ return b / (1.0 + dot(a, a)); });
      fv::cuBinary(gravity, centrifugal, centrifugal, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return fmin(10.0*a, fmax(-0.9*a, b)); });
    }

    //calculate advection
    fv::cuAdvectionNSWEs2ndRobustCurv(gravity, centrifugal, h, z, u, h_grad, eta_grad, u_grad, h_advection, hU_advection);

    //multiply advection with -1
    fv::cuUnaryOn(h_advection, [] __device__ (Scalar& a) -> Scalar{return -1.0*a;});
    fv::cuUnaryOn(hU_advection, [] __device__ (Vector& a) -> Vector{return -1.0*a;});

    fv::cuEulerIntegrator(h_old, h_advection_old, time_controller.dt() / 2.0, time_controller.current());
    fv::cuEulerIntegrator(hU_old, hU_advection_old, time_controller.dt() / 2.0, time_controller.current());
    fv::cuEulerIntegrator(h_old, h_advection, time_controller.dt() / 2.0, time_controller.current());
    fv::cuEulerIntegrator(hU_old, hU_advection, time_controller.dt() / 2.0, time_controller.current());

    //Write updated values
    fv::cuBinaryOn(h, h_old, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return b; });
    fv::cuBinaryOn(hU, hU_old, [] __device__(Vector& a, Vector& b) -> Vector{ return b; });

    hU.update_boundary_values();
    hU.update_time(time_controller.current(), time_controller.dt());

    fv::cuBinaryOn(hU, h, filter);
    fv::cuUnaryOn(h, filter1);

    //add centrifugal force to gravity to obtain gross vertical acceleration
    fv::cuBinary(gravity, centrifugal, gross_acc, [] __device__(Scalar& a, Scalar& b) -> Scalar{ return a + b; });

    //Calculating friction, use old h for depth, new hU for direction
    switch (rheology_type){
    case 1:
      fv::cuFrictionMCPlastic(time_controller.dt(), gross_acc, miu1, cohesion, rho, h, hU, z_grad, friction_force);
      break;
    case 2:
      fv::cuFrictionMuI(time_controller.dt(), gross_acc, miu1, miu2, _beta, _L, h, hU, z_grad, friction_force);
      break;
    case 3:
      fv::cuFrictionLucas(time_controller.dt(), gross_acc, miu1, miu2, _U, h, hU, z_grad, friction_force);
      break;
    }
    fv::cuEulerIntegrator(hU, friction_force, time_controller.dt(), time_controller.current());


    h.update_boundary_values();
    h.update_time(time_controller.current(), time_controller.dt()); 


    //forwarding the time
    time_controller.forward();
    time_controller.updateByCFL(gravity,h,hU);

    //print current time
    std::cout << time_controller.current() << std::endl;

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float elapsed_time = 0.0;
    cudaEventElapsedTime(&elapsed_time, start, stop);
    total_runtime += elapsed_time;

    if (time_controller.current() >= t_out){
      raster_writer.write(h, "h", t_out);
      raster_writer.write(eta, "elev", t_out);
      cuSimpleWriterLowPrecision(h, "h", t_out);
      cuSimpleWriterLowPrecision(hU, "hU", t_out);
      cuSimpleWriterLowPrecision(eta, "elev", t_out);
      t_out += dt_out;
    }

  } while (!time_controller.is_end());

  std::cout << "Total runtime " << total_runtime << "ms" << std::endl;

  return 0;

}

