#ifndef _1DGRID_H
#define _1DGRID_H

#include<vector>
#include<algorithm>
#include<iostream>
#include<cassert>

///
/// Create a 1D grid for advection
/// create a grid for the fluid in the domain from [xmin,xmax] including ghost points. Allow number of ghost points/cells, ng, on each end. The resolution of the domain is determined by nx, the numebr of points or cells. 
///

class _1DGrid{
  
private:
  ///
  ///physical xmin of the domain excluding ghost point
  ///
  
  double xmin; 

  ///
  ///physical xmax of the domain excluding ghost point
  ///
  
  double xmax; 

  int ng;  //number of ghost points on each end
  int nx;  //domain points excluding ghost points
  
  int ilo; //lower bound point of the domain excluding ghost point
  int ihi; //high bound point of the domain excluding ghost point

  double dx;  //separation between each point

  std::vector<double> state;  //current state, carry the state for each physical coordinate
  std::vector<double> state_init; //initial state
  std::vector<double> temp_state;
  std::vector<double> x; //physical x-coordinate including ghost points

  
public:
  _1DGrid(int _nx, int _ng= 2, double _xmin= 0.0, double _xmax= 1.0)
    :xmin{_xmin},
     xmax{_xmax}, 
     ng{_ng},
     nx{_nx},

     ilo{_ng},
     ihi{_ng+_nx-1},

     dx{(_xmax-_xmin)/(_nx-1)},
     state(_nx+2*_ng, 0.0),
     state_init(_nx+2*_ng, 0.0),
     temp_state(_nx+2*_ng, 0.0)
  {
    assert (_nx > 0);
    assert (_ng > 0);
    
    // initialize physical coordinate, centered cell 
    for (int i =0; i< _nx+2*_ng; ++i){
      x.push_back(_xmin + (i - _ng+0.5)*dx);
    }
  }

  // Returns a scratched array, 0 array
  std::vector<double> scratch_array();

  // invoke the boundary condition to update ghost cell
  void fill_BCs_diff();
  //  void fill_BCs_vol();

  // set initial condition
  void set_init(std::vector<double> init_vec);
    
  std::vector<double> get_state();
  
  // friend std::ostream& operator<< (std::ostream& os, const FDGrid& grid);

  // make advection_solver class a friend of FDGrid
   friend class advection_solver;
};

#endif
