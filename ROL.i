%module ROL
%{
#include "ROLVector.h"
#include "ROLObjective.h"
using namespace ROL;
%}

%include <ROL_StdVector.hpp>
%template(StdVectorDouble) ROL::StdVector<double>;
%include "ROLVector.h"

%include <ROL_Objective.hpp>
%template(ObjectiveDouble) ROL::Objective<double>;
%include "ROLObjective.h"
