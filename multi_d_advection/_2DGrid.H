#ifndef _2DGRID_H
#define _2DGRID_H

#include "_2DArray.H" 
#include <algorithm>

// Create a 2D grid for advection

class _2DGrid{
// create a grid for the fluid in the domain from [xmin,xmax] including ghost points. Allow number of ghost points/cells, ng, on each end. The resolution of the domain is determined by nx, the numebr of points or cells. 
  
private:
  double xmin; //physical xmin of the domain excluding ghost point
  double xmax; //physical xmax of the domain excluding ghost point
  
  double ymin;
  double ymax;

  int ng;  //number of ghost points on each end, assume the same in x, y

  int nx;  //domain points excluding ghost points
  int ny;
  
  int ilo_x; //lower bound point of the domain excluding ghost point
  int ihi_x; //high bound point of the domain excluding ghost point

  int ilo_y;
  int ihi_y;
  
  double dx;  //separation between each point
  double dy;
  
  _2DArray state;  //current state, carry the state for each physical coordinate
  _2DArray state_init; //initial state
  _2DArray temp_state;

  std::vector<double> x; //physical x-coordinate including ghost points
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
     state_init(_ny+2*_ng, _nx+2*_ng, 0.0),
     temp_state(_ny+2*_ng, _nx+2*_ng, 0.0)
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

  // Returns a scratched array, 0 array
  _2DArray scratch_array();

  // invoke the boundary condition to update ghost cell
  void fill_BCs_diff();
  //  void fill_BCs_vol();

  // set initial condition
  void set_init(_2DArray init_vec);

  inline _2DArray get_state(){
    return state;
  }

  inline std::vector<double> get_x(){
    return x;
  }

  inline std::vector<double> get_y(){
    return y;
  }
  // make advection_solver class a friend of FDGrid
   friend class advection_solver;
};

#endif
