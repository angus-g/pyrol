
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

  /// Implementation of pure virtual from c++
  double value(const ROL::Vector<double>& x, double& tol)
  {
    // FIXME: call to valuezz instead
    //    return valuezz(xxx, tol);
    return 0;
  }

  /// Python/SWIG implementation
  double valuezz(const ROL::MyVector& x, double& tol)
  {
    std::cout << "tol = " << tol << "\n";
    std::cout << "size(x) = " << x.dimension() << "\n";
    return 0.23254;
  }

  };

}
