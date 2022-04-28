#include"_2DGrid.H"
#include<iomanip>

_2DArray _2DGrid::scratch_array(){
  
  _2DArray v_out(nx + 2*ng, ny+2*ng, 0.0);
  return v_out;
}

// This seems to be better
void _2DGrid::fill_BCs_diff(bool x_update){
  if (x_update){
    for (int j = 0; j < ihi_y ; ++j){
      for (int n = 0; n < ng; ++n)
	{
	  state(j, ilo_x-n-1) = state(j, ihi_x-n-1);
	  state(j, ihi_x+n+1) = state(j, ilo_x+n+1);
	}
    }
  }
  else{
    for (int i = 0; i < ihi_x ; ++i){
      for (int n = 0; n < ng; ++n)
	{
	  state(ilo_y-n-1, i) = state(ihi_y-n-1, i);
	  state(ihi_y+n+1, i) = state(ilo_y+n+1, i);
	}
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

