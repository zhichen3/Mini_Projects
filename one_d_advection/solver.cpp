#include "solver.H"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>

// Define a helper function minmod:
double minmod(double a, double b){
  
  if ( (std::abs(a) < std::abs(b)) && (a * b > 0.0)){
    return a;
  }
  else if ( (std::abs(b) < std::abs(a)) && (a * b > 0.0)){
    return b;
    }
  else{
    return 0.0;
  }
}


void advection_solver::solve(){  
  // Write data to file
   std::ofstream of;
   std::string init_filename{"init_cond.dat"};
   std::string fin_filename{};
   fin_filename = method + "_" + slope_method + "_fin.dat";

   // Initialize initial state
  grid.set_init(init_cond(grid.x));
  t = 0.0;
  dt = dt_init;
  C = u*dt/grid.dx;
    
  std::cout << "Integration starting, writing to file...." << std::endl;

  write_file(init_filename);

  // Call the integration method
  if (method == "ftcs")
    {
      ftcs();
    }
  else if (method == "upwinding")
    {
      upwinding();
    }
  else if (method == "predictor_corrector")
    {
      predictor_corrector();
    }
  else if (method == "method_of_lines")
    {
      method_of_lines();
    }
  
  std::cout << "Integration ended, writing to file...." << std::endl;
  
  write_file(fin_filename);
  
  if (num_periods - std::trunc(num_periods) == 0.0){
    if (method == "ftcs" || method == "upwinding"){
      std::cout << "Relative Error of " << method << " Method is " << find_error() << "%" <<  std::endl;
    }
    else{
      std::cout << "Relative Error of " << method << " Method with " << slope_method << " slope is " <<  find_error() << "%" <<  std::endl;
    }
    
  }

}





void advection_solver::ftcs(){
// Integrate the system using ftcs method

  // create default zero array
  std::vector<double> state_new = grid.scratch_array();

  // Make sure boundary condition is fulfilled.
  grid.fill_BCs_diff();
  
  while (t < tmax) {
    //take care of last step where t+dt > tmax
    if ( t+dt > tmax){
      dt = tmax - t;
      C = u*dt/grid.dx;
    }
        
    // update states inside physical domain
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      state_new[i] = grid.state[i] - 0.5*C*(grid.state[i+1] - grid.state[i-1]);
    }
    
    grid.state = state_new;
    // update boundary ghost cells using updated physical domain 
    grid.fill_BCs_diff();
    t += dt;
  }
}






void advection_solver::upwinding(){
  // Integrate the system using upwinding method, the same as predictor_corrector with piecewise_const slope.

  // create default zero array
  std::vector<double> state_new = grid.scratch_array();
  // Make sure boundary condition is fulfilled.
  grid.fill_BCs_diff();
  
  while (t < tmax){
    //take care of last step where t+dt > tmax
    if ( t+dt > tmax){
      dt = tmax - t;
      C = u*dt/grid.dx;
    }

    // state update, loop over physical domain
    // for (int i = grid.ilo-grid.ng+1; i < grid.ihi+grid.ng ; ++i){ 
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      state_new[i] = grid.state[i] - C*(grid.state[i] - grid.state[i-1]);
    }
    
    grid.state = state_new;
    t += dt;
    grid.fill_BCs_diff();
  }
}





void advection_solver::predictor_corrector()
{

  std::vector<double> state_new = grid.scratch_array();
  grid.fill_BCs_diff();
  
  while (t < tmax){
    if (t+dt > tmax){
      dt = tmax - t;
      C = u*dt/grid.dx;
    }

    // for (int i = grid.ilo-grid.ng+1; i < grid.ihi+grid.ng; ++i){
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      state_new[i] = grid.state[i] + dt * rhs(i);
    }
    
    grid.state = state_new;
    grid.fill_BCs_diff();
    t += dt;
  }
}





void advection_solver:: method_of_lines()
{
  std::vector<double> k1 = grid.scratch_array();
  std::vector<double> k2 = grid.scratch_array();
  std::vector<double> k3 = grid.scratch_array();
  std::vector<double> k4 = grid.scratch_array();
  std::vector<double> init_state = grid.state;

  grid.fill_BCs_diff();
  
  while (t < tmax){
    if (t + dt > tmax){
      dt= tmax - t;
      C = u*dt/grid.dx;
    }
    
    // RK4

    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      k1[i] = rhs(i);
    }

    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      grid.state[i] = init_state[i] + 0.5*dt*k1[i];
    }
    
    grid.fill_BCs_diff();
    
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      k2[i] = rhs(i);
    }

    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      grid.state[i] = init_state[i] + 0.5*dt*k2[i];
    }

    grid.fill_BCs_diff();
    
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      k3[i] = rhs(i);
    }

    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      grid.state[i] = init_state[i] + dt*k3[i];
    }

    grid.fill_BCs_diff();
    
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      k4[i] = rhs(i);
    }
    
    for (int i = grid.ilo; i < grid.ihi+1; ++i){
      grid.state[i] = init_state[i] + dt/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    
    // RK2
     
    // for (int i = grid.ilo; i < grid.ihi+1; ++i){
    //   k1[i] = rhs(i);
    // }
    
    // for (int i = grid.ilo; i < grid.ihi+1; ++i){
    //   grid.state[i] = init_state[i] + 0.5*dt*k1[i];
    // }
    
    // grid.fill_BCs_diff();
    
    // for (int i = grid.ilo; i < grid.ihi+1; ++i){
    //   k2[i] = rhs(i);
    // }

    // for (int i = grid.ilo; i < grid.ihi+1; ++i){
    //   grid.state[i] = init_state[i] + dt*k2[i];
    // }

    grid.fill_BCs_diff();
    t += dt;
    init_state = grid.state;

  }
}




double advection_solver:: rhs(int ind){
  double rhs{};
  rhs = -C/dt*(riemann_selector(ind)-riemann_selector(ind-1));
  return rhs;
}






double advection_solver:: riemann_selector(int ind){
  // Returns the correct equation for the riemann problem:

  if (method == "predictor_corrector")
    {
      if (u > 0.0)
	{
	  double al{};
	  al = grid.state[ind] + 0.5*grid.dx*(1.0 - dt*u/grid.dx) * slope(ind);
	  return al;
	}
      else
	{
	  double ar{};
	  ar = grid.state[ind+1] - 0.5*grid.dx*(1.0 + dt*u/grid.dx) * slope(ind+1);
	  return ar;
	}
    }

  else if (method == "method_of_lines")
    
    {
      if (u > 0)
	{
	  double al{};
	  al = grid.state[ind] + 0.5*grid.dx*slope(ind);
	  return al;
	}
      else
	{
 	  double ar{};
	  ar = grid.state[ind+1] - 0.5*grid.dx*slope(ind+1);
	  return ar;
	}
    }
  else
    {
      return 1;
    }
}




double advection_solver::slope(int ind)
// Choose different slope mode
{
  double slope{};
  
  // Upwinding is the same as predictor_corrector with piecewise_const slope
  if (slope_method == "piecewise_const"){
    slope = 0.0;
  }
  
  else if (slope_method == "centered"){
    slope = 0.5 * (grid.state[ind+1] - grid.state[ind-1])/grid.dx;
  }

  else if (slope_method == "minmod"){
    double left_slope{};
    double right_slope{};
    left_slope = (grid.state[ind] - grid.state[ind-1])/grid.dx;
    right_slope = (grid.state[ind+1] - grid.state[ind])/grid.dx;
    slope = minmod(left_slope,right_slope);
  }

  else if (slope_method == "MC"){
    double left_slope{};
    double centered_slope{};
    double right_slope{};
    right_slope = 2.0*(grid.state[ind+1]-grid.state[ind])/grid.dx;
    centered_slope = 0.5*(grid.state[ind+1]-grid.state[ind-1])/grid.dx;
    left_slope = 2.0*(grid.state[ind] - grid.state[ind-1])/grid.dx;
    
    slope = minmod(centered_slope, minmod(left_slope,right_slope));
  }
    
  return slope;
}




 
void advection_solver:: print_state(){

  for (auto g: grid.state){
    std::cout << std::setw(15) << g;
  }
  std::cout << std::endl;
}




double advection_solver:: find_error(){
  // Find the relative error between the initial and final state

  //Make sure it evolved over integer period cycles.
  assert(num_periods - std::trunc(num_periods)==0);
  
  double err{0.0};
  double init{0.0};
  double rel_err{};
  
  for (int i = 0; i < static_cast<int>(grid.state.size()); ++i){
    err += std::pow(grid.state[i] - grid.state_init[i],2);
    init += std::pow(grid.state_init[i],2);
  }
  rel_err = std::sqrt(err/init)*100.0;
  
  return rel_err;
}



void advection_solver::write_file(std::string fname){
  std:: ofstream of;
  of.open(fname);
  for (int i = grid.ilo; i < grid.ihi+1; ++i){
    of << grid.x[i] << std::setw(15)
       << grid.state[i] << std::endl;
  }
  
  of.close();
}
