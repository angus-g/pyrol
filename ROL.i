%module ROL
%{
#include "ROLVector.h"
#include "ROLObjective.h"
using namespace ROL;
%}

%include <ROL_StdVector.hpp>
%template(StdVectorDouble) ROL::StdVector<double>;
%include "ROLVector.h"

// A typemap for ROL::Vector input
%typemap(in) ROL::Vector<double> const&
{
  std::cout << "Hello\n";
  auto ptr = $input;
  //   $1 = dynamic_cast<ROL::Vector<double>>($input);
}

// Apply 'double& tol'
%apply double& INPUT { double& tol };

%feature("notabstract") MyObjective;
%include <ROL_Objective.hpp>
%template(ObjectiveDouble) ROL::Objective<double>;
%include "ROLObjective.h"
