#ifndef SOLVER_H
#define SOLVER_H

#include "_2DGrid.H"
#include <functional>
#include <string>

///
/// Advection solver class
///
class Advection_solver{
private:
  ///
  /// Creating a 2D grid using _2DGrid class
  ///
  _2DGrid grid;

  ///
  /// velocity of advection in x-direction
  ///
  double u;

  ///
  /// velocity of advection in y-direction, since y axis counts from top to bottom, need to negate the input y velocity
  ///
  double v;

  ///
  /// Courant number
  ///
  double C;
  
  //  double num_periods;

  ///
  /// Current time of integration
  ///
  double t;

  ///
  /// Number of periods for advection. If whole number is given, error between the initial and final state can be calculated.
  ///
  double num_periods;
  
  ///
  /// function for initial condition that takes in physical domain x and y and returns the initial advection state
  ///
  std::function<_2DArray(const std::vector<double>&, const std::vector<double>&)> init_cond;

  ///
  /// Advection method of choice
  ///
  std::string method;

  ///
  /// slope method of choice
  ///
  std::string slope_method;
  
  ///
  /// Max time of integration
  ///
  double tmax;

  ///
  /// timestep
  ///
  double dt;

  ///
  /// initial timestep
  ///
  double dt_init;

  ///
  /// Strang-splitting method. Requires C<1 for stability
  ///
  void split();
  double rhs(const int j, const int i, const bool x_update);
  double riemann(const int j, const int i, const bool x_update);
  double slope(const int j, const int i, const bool x_update);
  
public:
  Advection_solver (_2DGrid _grid, double _u, double _v, double _C, double _num_periods, std::function<_2DArray(const std::vector<double>&, const std::vector<double>&)> _init_cond, const std::string& _method="split", const std::string& _slope_method="centered")
    : grid{_grid},
      u{_u},
      v{-_v},
      C{_C},
      t{0.0},
      num_periods{_num_periods},
      init_cond{_init_cond},
      method{_method},
      slope_method{_slope_method}
  {
    assert (_method == "split" || _method == "unsplit" || _method == "method_of_lines" );
    assert(_slope_method == "centered" || _slope_method == "piecewise_const" || _slope_method == "MC" || _slope_method == "minmod");
    
    assert(grid.ng > 1);

    double Tx = (grid.xmax - grid.xmin)/std::abs(_u);
    double Ty = (grid.ymax - grid.ymin)/std::abs(_v);

    if (Tx < Ty){
      tmax = _num_periods * Ty;
    }
    else{
      tmax = _num_periods * Tx;
    }

    
    dt = _C/(std::abs(_u)/grid.dx + std::abs(_v)/grid.dy); 
    dt_init = dt;
  }

  // general routine for solving 2D advection
  void solve();

  // print solution
  void print_state();

  
  // Find the error between the initial and final state for a square physical domain, same velocity for u and v, and integer number of periods.
  //  double find_error();

  // Write state data to file.
  void write_file(const std::string& fname);

};


#endif
