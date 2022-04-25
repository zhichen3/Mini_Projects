#ifndef SOLVER_H
#define SOLVER_H

#include "_2DGrid.H"
#include <functional>
#include <string>

class Advection_solver{
private:
  _2DGrid grid;
  double u;
  double v;
  double C;
  //  double num_periods;
  double t;
  double tmax;
  
  std::function<_2DArray(const std::vector<double>&, const std::vector<double>&)> init_cond;

  std::string method;
  std::string slope_method;

  double dt;
  double dt_init;

  // Requires C<1 for stability
  void split();

  
public:
  Advection_solver (_2DGrid _grid, double _u, double _v, double _C, double _tmax, std::function<_2DArray(const std::vector<double>&, const std::vector<double>&)> _init_cond, const std::string& _method="split", const std::string& _slope_method="centered")
    : grid{_grid},
      u{_u},
      C{_C},
      t{0.0},
      tmax{_tmax},
      init_cond{_init_cond},
      method{_method},
      slope_method{_slope_method}
  {
    assert (_method == "split" || _method == "unsplit" || _method == "method_of_lines" );
    assert(_slope_method == "centered" || _slope_method == "piecewise_const" || _slope_method == "MC" || _slope_method == "minmod");
    
    assert(grid.ng > 1);

    dt = _C/(_u/grid.dx + _v/grid.dy); 
    dta_init = dt;
  }

  // general routine for solving 2D advection
  void solve();

  // print solution
  void print_state();

  // Find the error between the initial and final state for a square physical domain, same velocity for u and v, and integer number of periods.
  double find_error();

  // Write state data to file.
  void write_file(std::string fname);
};


#endif
