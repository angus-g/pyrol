
#include <ROL_StdVector.hpp>
#include <vector>

namespace ROL
{
  class MyVector : public StdVector<double>
  {
  public:
  MyVector(int n) : StdVector<double>(Teuchos::rcp(new std::vector<double>(n)))
    {
    }

  };

}
