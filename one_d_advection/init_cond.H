#ifndef INIT_COND_H
#define INIT_COND_H

#include<vector>
#include<cmath>

// This file stores some simple initial conditions 

// tophat initial condition
// vec is the physical domain
std::vector<double> tophat(const std::vector<double>& domain){
  std::vector<double> state(domain.size(), 0.0);
  for (int i = 0; i < static_cast<int>(domain.size()); ++i){
    if (domain[i] > 1.0/3.0 && domain[i] < 2.0/3.0){
      state[i] = 1.0;
    }
  }

  return state;
}

std::vector<double> sine(const std::vector<double>& domain){
  std::vector<double> state(domain.size(), 1.0);
  const double pi = 3.141592654;

  for (int i = 0; i < static_cast<int>(domain.size()); ++i){
    state[i] += 0.5 * std::sin(2.0*pi*domain[i]);
  }
  
  return state;
}


std::vector<double> gaussian(const std::vector<double> &domain){
  std::vector<double> state(domain.size(), 0.0);
  const double sigma = 0.15;
  for (int i = 0; i < static_cast<int>(domain.size()); ++i){
    state[i] = std::exp(-0.5*std::pow(domain[i]-0.5,2)/std::pow(sigma,2));
  }
  return state;
}


#endif
