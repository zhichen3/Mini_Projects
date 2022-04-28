#include "init_cond.H"
#include "_2DGrid.H"
#include "solver.H"
#include <iostream>
#include <fstream>

int main() {

  int nx{64};
  int ny{64};
  int ng{2};
  _2DGrid grid(nx,ny,ng);
  _2DArray tophat_state = tophat(grid.get_x(), grid.get_y());
  _2DArray gaussian_state = gaussian(grid.get_x(), grid.get_y());

  // std::ofstream of;
  // of.open("gaussian.dat");
  // for (int j = ng; j < ng+ny; ++j){
  //   for (int i=ng; i < ng+nx; ++i){
  //     of << gaussian_state(j,i) << " ";
  //   }
  //   of << std::endl;
  // }
  // of.close();


  // of.open("tophat.dat");
  // for (int row =0; row < tophat_state.nrows(); ++row){
  //   for (int col=0; col < tophat_state.ncols(); ++col){
  //     of << tophat_state(row,col) << " ";
  //   }
  //   of << std::endl;
  // }
  // of.close();


  
  double u{1.0};
  double v{0.5};
  double C{0.3};
  double num_periods{0.32};

  Advection_solver split_centered(grid, u, v, C, num_periods, gaussian, "split", "centered");
  split_centered.solve();
  
  
}
