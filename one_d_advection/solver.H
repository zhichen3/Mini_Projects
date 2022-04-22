#ifndef SOLVER_H
#define SOLVER_H

#include "_1DGrid.H"
#include <functional>
#include <cmath>
#include <string>

class advection_solver{

private:
  _1DGrid grid;
  double u;
  double C;
  double num_periods;
  double t;
  
  std::function<std::vector<double>(std::vector<double>&)> init_cond;

  std::string method;
  std::string slope_method;
  
  double dt;
  double dt_init;
  double tmax;

  // ftcs method
  void ftcs();
  
  // upwinding method, requires C<=1 for stability
  void upwinding();

  // predictor_corrector method:
  // if choose slope to be piecewise_const, then its upwinding method
  void predictor_corrector();

  // method of lines method,
  void method_of_lines();

  // RHS of equation da/dt
  double rhs(int ind);
  
  // Riemann Problem Selector
  double riemann_selector(int ind);

  // Find slope between cells using center-cell
  double slope(int ind);

public:
  advection_solver (_1DGrid _grid, double _u, double _C, std::function<std::vector<double>(const std::vector<double>& )> _init_cond, double _num_periods=1.0, const std::string& _method="upwinding", const std::string& _slope_method="centered")
    : grid{_grid},
      u{_u},
      C{_C},
      num_periods{_num_periods},
      t{0.0},  
      init_cond{_init_cond},
      method{_method},
      slope_method{_slope_method}
  {
    // assert correct advection methods
    assert(_method == "ftcs" || _method == "upwinding" || _method == "predictor_corrector" || _method == "method_of_lines");

    // assert correct slope methods
    assert(_slope_method == "centered" || _slope_method == "piecewise_const" || _slope_method == "MC" || _slope_method == "minmod");
    
    //make sure at least 2 ghost cells for second order methods
    if (_method == "predictor_corrector" || _method == "method_of_lines")
      {
	assert(grid.ng > 1);
      }
    
    dt = _C*grid.dx/_u;
    dt_init = dt;
    tmax = _num_periods*(grid.xmax-grid.xmin)/std::abs(_u);
  }

  // general routine for solving advection equation.
  void solve();
     
  // print solution
  void print_state();

  // Find the error between initial and final state
  double find_error();

  // write state data to file 
  void write_file(std::string fname);
};


#endif
