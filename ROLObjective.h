
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

  /// Implementation of pure virtual from c++ (ignore in SWIG)
  double value(const ROL::Vector<double>& x, double& tol)
  {
    return value(dynamic_cast<const ROL::MyVector&>(x), tol);
  }

  /// Python/SWIG implementation
  double value(const ROL::MyVector& x, double& tol)
  {
    std::cout << "tol = " << tol << "\n";
    std::cout << "size(x) = " << x.dimension() << "\n";
    return 0.23254;
  }

  };

}
