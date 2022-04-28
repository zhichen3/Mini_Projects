#include "solver.H"
#include "_1DGrid.H"
#include "init_cond.H"

int main()
{
  int nx{64};
  int ng{2};
  double xmin{0.0};
  double xmax{1.0};
  _1DGrid g(nx,ng,xmin,xmax);

  double u{1.0};
  double C{0.5};
  double num_periods{3.0};
  
  advection_solver tophat_pc_pc(g,u,C, tophat, num_periods, "predictor_corrector", "MC", true);
  tophat_pc_pc.solve();
}
