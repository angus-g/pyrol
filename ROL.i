%module ROL
%{
#include "ROLVector.h"
#include "ROLObjective.h"
using namespace ROL;
%}

// -- ROL::Vector<double> --

%include <ROL_StdVector.hpp>
%template(StdVectorDouble) ROL::StdVector<double>;
%include "ROLVector.h"

// -- ROL::Objective<double> --

// Apply 'double& tol'
%apply double& INPUT { double& tol };

%feature("notabstract") MyObjective;
%include <ROL_Objective.hpp>
%template(ObjectiveDouble) ROL::Objective<double>;
%ignore ObjectiveDouble::value(const ROL::Vector<double>& x, double& tol);

%include "ROLObjective.h"

// -- ROL::Algorithm<double> --
 //
 // TODO
