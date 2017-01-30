%module ROL
%{
#include "ROLVector.h"
#include "ROLObjective.h"
using namespace ROL;
%}

%include <ROL_StdVector.hpp>
%template(StdVectorDouble) ROL::StdVector<double>;
%include "ROLVector.h"

// // A typemap for ROL::Vector input
// %typemap(in) ROL::Vector<double> const&
// {
//   std::cout << "Hello\n";

//  //  $input = $1;


//   //  std::cout << "$input = " << $input << "\n";
//   //  $1 = reinterpret_cast< ROL::MyVector * >($input);
//   //  std::cout << "$input = " << $1 << "\n";
//   auto ptr = std::make_shared<ROL::MyVector>(200);
//   $1 = ptr.get();
// }

// Apply 'double& tol'
%apply double& INPUT { double& tol };

%feature("notabstract") MyObjective;
%include <ROL_Objective.hpp>
%template(ObjectiveDouble) ROL::Objective<double>;
%include "ROLObjective.h"
