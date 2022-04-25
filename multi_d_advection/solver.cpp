#include "solver.H"

void Advection_solver::solver(){
  std::ofstream of;
  std::string init_filename{"init_cond.dat"};
  std::string fin_filename{};
  fin_filename = method + "_" + slope_method + "_fin.dat";

  //set initial condition
  grid.set_init(init_cond(grid.x,grid.y));
  t = 0.0;
  dt = dt_init;
  C = dt*(u/grid.dx + v/grid.dy);

  std::cout << "Integration starting, writing to file ..." << std::endl;

  write_file(init_filename);

  if (method == "split")
    {
      split();
    }
  else if (method == "unsplit")
    {
      unsplit();
    }
  else if (method == "method_of_lines")
    {
      method_of_lines();
    }

  std::cout << "Integration ended, writing to file ..." << std::endl;

  write_file(fin_filename);
}



void Advection_solver::write_file(std::string fname){
  std::ofstream of;
  of.open(fname);

  for (int i = grid.ilo_y; j < grid.ilo_y+1; ++i){
    for (int j = grid.ilo_x; i < grid.ilo_x+1; ++j){
      of << state(i,j) << " ";
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
    
    // Update state_new:    
    if x_update{

	for (int i = grid.ilo_y; j < grid.ihi_y+1; ++i){
	  for (int j = grid.ilo_x; i < grid.ihi_x+1; ++j){
	    state_temp(i,j) = grid.state(i,j) - dt*u/grid.dx* rhs(grid.state(i,j));
	  }
	}
	
	for (int i = grid.ilo_y; j < grid.ihi_y+1; ++j){
	  for (int j = grid.ilo_x; i < grid.ihi_x+1; ++i){
	    grid.state(i,j) = state_temp(i,j) + dt*v/grid.dy*rhs(temp_state(i,j));
	  }
	}	
      }
    
    else {
      
      for (int i = grid.ilo_y; j < grid.ihi_y+1; ++i){
	for (int j = grid.ilo_x; i < grid.ihi_x+1; ++j){
	  state_temp(i,j) = grid.state(i,j) - dt*v/grid.dy* rhs(grid.state(i,j));
	}
      }
      
      for (int i = grid.ilo_y; j < grid.ihi_y+1; ++i){
	for (int j = grid.ilo_x; i < grid.ihi_x+1; ++j){
	  grid.state(i,j) = state_temp(i,j) + dt*u/grid.dx*rhs(temp_state(i,j));  
	}
      }
      
    }
    
    
    grid.state = state_new;
    grid.fill_BCs_diff();
    t += dt;
    x_update = !x_update; 
  }
}

