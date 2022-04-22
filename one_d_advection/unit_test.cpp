#include"solver.H"
#include"_1DGrid.H"
#include"init_cond.H"
#include<algorithm>

int main(){ 
  int nx{128};
  int ng{2};
  double xmin{0.0};
  double xmax{1.0};
  _1DGrid g(nx,ng,xmin,xmax);

  ///
  /// Testing _1DGrid
  ///

  // Test _1DGrid.cratch_array()
  std::vector<double> vec(nx+2*ng, 0.0);
  for (int i = 0; i < nx+2*ng; ++i){
    assert(vec[i] == g.scratch_array()[i]);
  }

  // Test fill_BCs_diff and set_init functions:
  double k{0.0};
  for (auto it = vec.begin(); it < vec.end(); ++it){
    *it += k;
    ++k;
  }
  g.set_init(vec);
  for (int i = 0; i < nx+2*ng; ++i){
    assert(g.get_state()[i] == vec[i]);
  }

  g.fill_BCs_diff();
  for (int n = 0; n < ng; ++n){
    assert(g.get_state()[ng-n-1] == g.get_state()[ng+nx-n-2]);
    assert(g.get_state()[ng+nx+n] == g.get_state()[ng+n+1]);
  }

  ///
  /// Test different advection solver
  ///
  double u{1.0};
  double C{0.5};
  double num_periods{3.0};

  ///
  /// Test linear advection solver with tophat initial cond
  ///

  // use lower value of C for ftcs method.
  advection_solver tophat_ftcs(g, u, 0.001, tophat, num_periods, "ftcs");
  tophat_ftcs.solve();  
  assert(tophat_ftcs.find_error() < 50.0);

  advection_solver tophat_upwinding(g, u, C, tophat, num_periods, "upwinding");
  tophat_upwinding.solve();  
  assert(tophat_upwinding.find_error() < 50.0);

  ///
  /// Test nonlinear advection solver with tophat inital cond
  ///

  // Predictor_protector method, for predictor_corrector with piecewise_cosnt slope should be the same as upwinding method
  advection_solver tophat_pc_pc(g, u, C, tophat, num_periods, "predictor_corrector", "piecewise_const");
  tophat_pc_pc.solve();  
  assert(tophat_pc_pc.find_error() -tophat_upwinding.find_error() < 0.001);

  advection_solver tophat_pc_centered(g, u, C, tophat, num_periods, "predictor_corrector", "centered");
  tophat_pc_centered.solve();  
  assert(tophat_pc_centered.find_error() < 50.0);  

  advection_solver tophat_pc_minmod(g, u, C, tophat, num_periods, "predictor_corrector", "minmod");
  tophat_pc_minmod.solve();  
  assert(tophat_pc_minmod.find_error() < 50.0);

  advection_solver tophat_pc_MC(g, u, C, tophat, num_periods, "predictor_corrector", "MC");
  tophat_pc_MC.solve();  
  assert(tophat_pc_MC.find_error() < 50.0);


  // Method of lines method:
  
  advection_solver tophat_ml_pc(g, u, C, tophat, num_periods, "method_of_lines", "piecewise_const");
  tophat_ml_pc.solve();  
  assert(tophat_ml_pc.find_error() < 50.0);

  advection_solver tophat_ml_centered(g, u, C, tophat, num_periods, "method_of_lines", "centered");
  tophat_ml_centered.solve();  
  assert(tophat_ml_centered.find_error() < 50.0);  

  advection_solver tophat_ml_minmod(g, u, C, tophat, num_periods, "method_of_lines", "minmod");
  tophat_ml_minmod.solve();  
  assert(tophat_ml_minmod.find_error() < 50.0);

  advection_solver tophat_ml_MC(g, u, C, tophat, num_periods, "method_of_lines", "MC");
  tophat_ml_MC.solve();  
  assert(tophat_ml_MC.find_error() < 50.0);

}

