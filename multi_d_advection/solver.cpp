#include "solver.H"
#include <iostream>
#include <fstream>

void Advection_solver::solve(){
 
  std::string init_filename{"init_cond.dat"};
  std::string fin_filename{};
  fin_filename = method + "_" + slope_method + "_fin.dat";

  //set initial condition
  grid.set_init(init_cond(grid.x, grid.y));
  t = 0.0;
  dt = dt_init;
  C = dt*(u/grid.dx + v/grid.dy);

  std::cout << "Integration starting, writing to file ..." << std::endl;

  write_file(init_filename);
  
  if (method == "split")
    {
      split();
    }
  // else if (method == "unsplit")
  //   {
  //     unsplit();
  //   }
  // else if (method == "method_of_lines")
  //   {
  //     method_of_lines();
  //   }


  std::cout << "Integration ended, writing to file ..." << std::endl;
   write_file(fin_filename);
   
}


void Advection_solver::print_state(){
  std::cout << grid.get_state() << std::endl;
}


void Advection_solver::write_file(const std::string& fname){
  std::ofstream of;
  of.open(fname);

  for (int j = grid.ilo_y; j < grid.ihi_y+1; ++j){
    for (int i = grid.ilo_x; i < grid.ihi_x+1; ++i){
      of << grid.state(j,i) << " ";
    }
    of << std::endl;
  }
  
  of.close();
}



void Advection_solver::split(){
 
  bool x_update = true; 
  _2DArray temp_state = grid.scratch_array();
  grid.fill_BCs_diff();

  
  while (t < tmax){
    if (t + dt > tmax){
      dt = tmax - t;
      C = dt*(u/grid.dx+v/grid.dy);
    }
    
    for (int k = 0; k < 2; ++k){
      for (int j = grid.ilo_y; j < grid.ihi_y+1; ++j){
	for (int i = grid.ilo_x; i < grid.ihi_x+1; ++i){
	  temp_state(j,i) = grid.state(j,i) + dt* rhs(j, i, x_update);
	}
      }
      
      grid.state = temp_state;
      grid.fill_BCs_diff();
      x_update = !x_update;
    }
    
    x_update = !x_update;
    t += dt;
   
  }
}

double Advection_solver:: rhs(const int j, const int i, const bool x_update){
  double rhs{};
  
  if (x_update){
       rhs = -u/grid.dx*(riemann(j, i, x_update) - riemann(j, i-1, x_update) );
  }
  else{
    rhs = -v/grid.dy *(riemann(j, i, x_update) - riemann(j-1, i, x_update) );
    
  }
  
  return rhs;
}


double Advection_solver:: riemann(const int j, const int i, const bool x_update){

  double res{};
  
  if (x_update){
    if (u > 0){
      res = grid.state(j,i) + 0.5*(grid.dx - u*dt)*slope(j, i, x_update);
    }

    else {
      res = grid.state(j,i+1) - 0.5*(grid.dx + u*dt)*slope(j, i+1, x_update);
    }

  }
  else{

    if (v > 0){
      res = grid.state(j,i) + 0.5*(grid.dy - v*dt)*slope(j, i, x_update);
    }

    else{
      res = grid.state(j+1,i) - 0.5*(grid.dy + v*dt)*slope(j+1, i, x_update);

    }

  }

  return res;
}


double Advection_solver::slope(const int j, const int i, const bool x_update){
  double slope{};
  
  if (slope_method == "piecewise_const"){
    slope = 0.0;
  }

  if (x_update){
    if (slope_method == "centered"){
      slope = 0.5*(grid.state(j,i+1) - grid.state(j,i-1))/grid.dx;
    }
  }
  else{
    if (slope_method == "centered"){
      slope = 0.5*(grid.state(j+1,i) - grid.state(j-1,i))/grid.dy;
    }
  }
  return slope;
}

			       
