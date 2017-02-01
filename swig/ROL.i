%module(directors="1") ROL
%{
#include "ROLVector.h"
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>
using namespace ROL;
%}

%include "exception.i"
%include "std_string.i"

// Typemap to convert a sequence to a Teuchos pointer to std::vector<double>
%typemap(in) Teuchos::RCP<std::vector<double>> const&
(Teuchos::RCP<std::vector<double>> tmp_ptr)
{

  if (!PySequence_Check($input))
  {
    SWIG_exception(SWIG_TypeError, "expected a sequence for argument $argnum");
  }

  // Get sequence length
  Py_ssize_t pyseq_length = PySequence_Size($input);
  tmp_ptr = Teuchos::rcp<std::vector<double>>(new std::vector<double>(pyseq_length));
  std::vector<double>& tmp_array = *tmp_ptr;

  for (int i = 0; i < pyseq_length; i++)
  {
    double value;
    PyObject *item = PySequence_ITEM($input, i);

    if(!SWIG_IsOK(SWIG_AsVal(double)(item, &value)))
    {
      Py_DECREF(item);
      SWIG_exception(SWIG_TypeError, "expected items of sequence to be of type " \
      		     "\"double\" in argument $argnum");
    }
    tmp_array[i] = value;
    Py_DECREF(item);
  }

  $1 = &tmp_ptr;
}

// -- ROL::Vector<double> --

%include <ROL_Vector.hpp>
%template(VectorDouble) ROL::Vector<double>;

%include <ROL_StdVector.hpp>
%template(StdVectorDouble) ROL::StdVector<double>;

%include "ROLVector.h"

// -- ROL::Objective<double> --

// Apply 'double& tol'
%apply double& INPUT { double& tol };

%feature("director") ROL::Objective<double>;

%include <ROL_Objective.hpp>
%template(ObjectiveDouble) ROL::Objective<double>;

// -- ROL::Algorithm<double> --

%typemap(in) Teuchos::ParameterList&
(Teuchos::ParameterList parlist)
{
  // Need to create Parameterlist from python input, somehow.

  parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Newton-Krylov");
  parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
  parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
  parlist.sublist("Status Test").set("Iteration Limit",100);

  $1 = &parlist;
}

%typemap(typecheck) Teuchos::ParameterList&
{
  // Always pass as true... FIXME
  $1 = 1;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) ROL::Vector< double > &
{
  // Always pass as true... FIXME
  $1 = 1;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) ROL::Objective<double>&
{
  // Always pass as true... FIXME
  $1 = 1;
}

%include <ROL_Algorithm.hpp>
%template(AlgorithmDouble) ROL::Algorithm<double>;
