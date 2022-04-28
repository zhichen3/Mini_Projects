#ifndef _2DGRID_H
#define _2DGRID_H

#include "_2DArray.H" 
#include <algorithm>

///
/// create a grid for the fluid in the domain from [xmin,xmax] including ghost points. Allow number of ghost points/cells, ng, on each end. The resolution of the domain is determined by nx, the numebr of points or cells. 
///
class _2DGrid{
  
private:

  ///
  /// Physical xmin of the domain excluding ghost point
  ///
  double xmin;

  ///
  /// Physical xmax of the domain excluding ghost point
  ///
  double xmax;

  ///
  /// Physical ymin of the domain excluding ghost point
  ///
  double ymin;

  ///
  /// Physical ymax of the domain excluding ghost point
  ///
  double ymax;

  ///
  /// Number of ghost points on each end, (assumes number of ghost points are the same in both x and y)
  ///
  int ng;

  ///
  /// Number of domain points in x direction excluding ghost points
  ///
  int nx;

  ///
  /// Number of domain points in y direction exlcuding ghost points
  ///
  int ny;

  ///
  /// Lower bound point of the domain, x-indice, excluding ghost points
  ///
  int ilo_x;

  ///
  /// high bound point of the domain, x-indice, excluding ghost point
  ///
  int ihi_x; 

  ///
  /// low bound point of the domain, y-indice, excluding ghost points
  ///
  int ilo_y;

  ///
  /// high bound point of the domain, y-indice,  excluding ghost point
  ///
  int ihi_y;

  ///
  /// Separation between each physical domain in the x-direction
  ///
  double dx;

  ///
  /// Separation between each physical domain in the y-direction
  ///
  double dy;

  ///
  /// Current state, carry the state for each physical coordinate
  ///
  _2DArray state;

  ///
  /// Initial state of advection
  _2DArray state_init; 
  //  _2DArray temp_state;

  ///
  /// Physical x-coordinate including ghost points
  ///
  std::vector<double> x;

  ///
  /// Physical y-coordinate including ghost points
  ///
  std::vector<double> y; 
  
public:
  _2DGrid(int _nx, int _ny, int _ng= 2, double _xmin= 0.0, double _xmax= 1.0, double _ymin=0.0, double _ymax=1.0)
    :xmin{_xmin},
     xmax{_xmax},
     ymin{_ymin},
     ymax{_ymax},
     
     ng{_ng},
     nx{_nx},
     ny{_ny},
     
     ilo_x{_ng},
     ihi_x{_ng+_nx-1},
     ilo_y{_ng},
     ihi_y{_ng+_ny-1},
     
     dx{(_xmax-_xmin)/(_nx-1)},
     dy{(_ymax-_ymin)/(_ny-1)},

     // set initial state to 0
     state(_ny+2*_ng, _nx+2*_ng, 0.0),
     state_init(_ny+2*_ng, _nx+2*_ng, 0.0)
     //     temp_state(_ny+2*_ng, _nx+2*_ng, 0.0)
  {
    assert (_nx > 0);
    assert (_ny > 0);
    assert (_ng > 0);
    
    // initialize physical coordinate, centered cell 
    for (int i =0; i< _nx+2*_ng; ++i){
      x.push_back(_xmin + (i - _ng+0.5)*dx);
    }
    
    for (int i =0; i< _ny+2*_ng; ++i){
      y.push_back(_ymin + (i - _ng+0.5)*dy);
    }
    
  }

  ///
  /// Returns a scratched 2D array, initialized to 0
  ///
  _2DArray scratch_array();

  ///
  /// invoke the boundary condition to update ghost cell
  ///
  void fill_BCs_diff(bool x_update);
  

  ///
  /// set initial condition
  ///
  void set_init(_2DArray init_vec);

  ///
  /// Returns the current advection state
  ///
  inline _2DArray get_state(){
    return state;
  }
  
  inline std::vector<double> get_x(){
    return x;
  }

  inline std::vector<double> get_y(){
    return y;
  }

  
  /// make advection_solver class a friend of 2DGrid
   friend class Advection_solver;
};

#endif
