#include "init_cond.H"
#include "_2DGrid.H"
#include <iostream>
#include <fstream>

int main() {

  int nx{64};
  int ny{64};
  _2DGrid grid(nx,ny);
  _2DArray tophat_state = tophat(grid.get_x(), grid.get_y());
  _2DArray gaussian_state = gaussian(grid.get_x(), grid.get_y());
  // std::cout << gaussian_state <<std::endl;
  // std::cout << std::endl;
  // std::cout << grid.get_state() << std::endl;

  std::ofstream of;
  of.open("gaussian.dat");
  for (std::size_t row =0; row < tophat_state.nrows(); ++row){
    for (std::size_t col=0; col < tophat_state.ncols(); ++col){
      of << gaussian_state(row,col) << " ";
    }
    of << std::endl;
  }
  of.close();


  of.open("tophat.dat");
  for (std::size_t row =0; row < tophat_state.nrows(); ++row){
    for (std::size_t col=0; col < tophat_state.ncols(); ++col){
      of << tophat_state(row,col) << " ";
    }
    of << std::endl;
  }
  of.close();

}
