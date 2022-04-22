#include"_1DGrid.H"
#include<iomanip>

std::vector<double> _1DGrid::scratch_array(){
  
  std::vector<double> v_out(nx + 2*ng, 0.0);
  return v_out;
}

// This seems to be better
void _1DGrid::fill_BCs_diff(){

  for (int n = 0; n < ng; ++n)
    {
      state[ilo-n-1] = state[ihi-n-1];
      state[ihi+n+1] = state[ilo+n+1];
    }
}


// void _1DGrid::fill_BCs_vol(){
//   for (int n = 0; n < ng; ++n)
//     {
//       state[ilo-n-1] = state[ihi-n];
//       state[ihi+n+1] = state[ilo+n];
//     }
// }


void _1DGrid::set_init(std::vector<double> init_vec){

  assert(init_vec.size()==state.size());
 
  state = init_vec;
  state_init = init_vec;
}

std::vector<double> _1DGrid::get_state(){
  return state;
}

// std::ostream& operator<< (std::ostream& os, const FDGrid& grid){
  
//   for (auto g : grid.x){
//     os << std::setw(12) << g;
//   }
//   os << std::endl;
//   return os;
//}
