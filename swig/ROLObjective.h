
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
  virtual double value(const ROL::Vector<double>& x, double& tol)
  {
    std::cout << "Error: please overload value\n";
    exit(-1);
    return 0.0; //value(dynamic_cast<const ROL::MyVector&>(x), tol);
  }

  // /// Python/SWIG implementation
  // double value(const ROL::MyVector& x, double& tol)
  // {
  //   std::cout << "tol = " << tol << "\n";
  //   std::cout << "size(x) = " << x.dimension() << "\n";
  //   return x.norm();
  // }

  };

}
