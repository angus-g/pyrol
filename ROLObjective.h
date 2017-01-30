
#include <ROL_Objective.hpp>
#include <vector>

namespace ROL
{
  class MyObjective : public Objective<double>
  {
  public:

  MyObjective()
    {
    }

  double value(const ROL::Vector<double>& x, double& tol)
  {
    return 0.23254;
  }

  };

}
