#ifndef INIT_COND_H
#define INIT_COND_H

#include "_2DArray.H"
#include <cmath>



_2DArray tophat(const std::vector<double>& x, const std::vector<double>& y){
  _2DArray state(static_cast<std::size_t>(y.size()), static_cast<std::size_t>(x.size()), 0.0);

  for (int j = 0; j < static_cast<int>(y.size()); ++j){
    for (int i = 0; i < static_cast<int>(x.size()); ++i){
      if ( std::pow(x[i]-0.5, 2) + std::pow(y[j]-0.5,2) < std::pow(1.0/3.0,2) ) {
	state(j,i) = 1.0;
      }
    }
  }
  return state;
}





_2DArray gaussian(const std::vector<double>& x, const std::vector<double>& y){
  _2DArray state(static_cast<int>(y.size()), static_cast<int>(x.size()), 0.0);
  const double sigma_x{0.15};
  const double sigma_y{0.15};
  
  for (int j = 0; j < static_cast<int>(y.size()); ++j){
    for(int i = 0; i < static_cast<int>(x.size()); ++i){
      state(j,i) = std::exp(-0.5*( std::pow(x[i]-0.5,2)/(sigma_x*sigma_x)+std::pow(y[j]-0.5,2)/(sigma_y*sigma_y) ));
    }
  }

  return state;
}




#endif
