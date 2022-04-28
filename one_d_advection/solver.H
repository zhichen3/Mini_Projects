#ifndef SOLVER_H
#define SOLVER_H

#include "_1DGrid.H"
#include <functional>
#include <cmath>
#include <string>

///
/// 1-D advection solver. Contains both linear and nonlinear advection solver.
///
class advection_solver{

private:
  ///
  /// 1DGrid class using _1DGrid.H
  ///
  _1DGrid grid;

  ///
  /// velocity of advection.
  ///Positive travels to the right. Negative travels to the left
  ///
  double u;

  ///
  /// Courant–Friedrichs–Lewy number, which is a a measure of what fraction of a zone we cross per timestep
  ///
  double C;

  ///
  /// Number of periods for advection
  ///
  double num_periods;

  ///
  /// current time for advection
  ///
  double t;

  ///
  /// A function that takes the physical domain, e.g. x, as argument that returns the initial state for advection.
  ///
  std::function< std::vector<double> (std::vector<double>&) > init_cond;

  ///
  /// Advection method of choice.
  /// For linear advection, choose between ftcs or upwinding
  /// For nonlinear advection, choose between predictor_corrector or method_of_lines
  ///
  std::string method;

  ///
  /// Slope method of choice for nonlinear advection methods.
  /// Choose between piecewise_const, centered, minmod, MC
  ///
  std::string slope_method;

  ///
  /// Set this to true to allow runtime advection states, which then can be used for doing runtime visualization.
  ///
  bool runtime_dat;
  
  ///
  /// current time step
  ///
  double dt;

  ///
  /// initial timestep
  ///
  double dt_init;

  ///
  /// end time of advection
  ///
  double tmax;
  
  ///
  /// Linear advection method: ftcs method. Smaller C works better
  ///
  void ftcs();

  ///
  /// Linear upwinding method, requires C<=1 for stability. If C=1 we have an exact solution
  ///
  void upwinding();

  ///
  /// Nonlinear advection method: predictor_corrector.
  /// if choose slope to be piecewise_const, then its the same as upwinding method.
  ///
  void predictor_corrector();

  ///
  /// Nonlinear advection method: method of lines method, uses RK4 integration.
  ///
  void method_of_lines();

  ///
  /// Returns the RHS of equation da/dt
  ///
  double rhs(int ind);

  ///
  /// Returns the result of Riemann Problem depending on which advection directio, i.e. depending on the sign of u.
  ///
  double riemann_selector(int ind);

  ///
  /// Returns the slope depending on the slope method of choice.
  ///
  double slope(int ind);

  
public:
  ///
  /// Constructing the advection solver
  ///
  advection_solver (_1DGrid _grid, double _u, double _C, std::function<std::vector<double>(const std::vector<double>& )> _init_cond, double _num_periods=1.0, const std::string& _method="upwinding", const std::string& _slope_method="centered",const bool _runtime_dat = false)
    : grid{_grid},
      u{_u},
      C{_C},
      num_periods{_num_periods},
      t{0.0},  
      init_cond{_init_cond},
      method{_method},
      slope_method{_slope_method},
      runtime_dat{_runtime_dat}
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
  
  ///
  /// It solves the advection equation with given advection method and slope if chose a nonlinear method.
  ///
  void solve();

  ///
  /// It prints the current state(solution)
  ///
  void print_state();

  ///
  /// Return the relative error between initial and final state
  /// It can be called only when num_of_periods is a whole number, i.e. making a full cycle of advection.
  /// It is used to test the validity of advection method, which should return back to the original state.
  ///
  double find_error();

  ///
  /// A procedure to  write state data to file
  ///
  void write_file(std::string fname);
};


#endif
