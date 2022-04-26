#include"_2DGrid.H"
#include<iomanip>

_2DArray _2DGrid::scratch_array(){
  
  _2DArray v_out(nx + 2*ng, ny+2*ng, 0.0);
  return v_out;
}

// This seems to be better
void _2DGrid::fill_BCs_diff(){
  for (int i = 0; i < ihi_y ; ++i){
    for (int n = 0; n < ng; ++n)
      {
	state(i, ilo_x-n-1) = state(i, ihi_x-n-1);
	state(i, ihi_x+n+1) = state(i, ilo_x+n+1);
      }
  }

  for (int j = 0; j < ihi_x ; ++j){
    for (int n = 0; n < ng; ++n)
      {
	state(ilo_y-n-1, j) = state(ihi_y-n-1, j);
	state(ihi_y+n+1, j) = state(ilo_y+n+1, j);
      }
  }
}

// void _1DGrid::fill_BCs_vol(){
//   for (int n = 0; n < ng; ++n)
//     {
//       state[ilo-n-1] = state[ihi-n];
//       state[ihi+n+1] = state[ilo+n];
//     }
// }


void _2DGrid::set_init(_2DArray init_vec){

  assert(init_vec.size() == state.size());
 
  state = init_vec;
  state_init = init_vec;
}

