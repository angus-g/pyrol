from hs.hs4 import run_HS4
from hs.hs13 import run_HS13
from hs.hs28 import run_HS28
from hs.hs29 import run_HS29
import ROL
from ROL.numpy_vector import NumpyVector
params_dict = {
    "Step": {
        "Moreau-Yosida Penalty":  {
            "Initial Penalty Parameter": 1.0,
            "Subproblem": {
                "Iteration Limit": 200,
            }
        },
        "Trust Region":
        {
            "Subproblem Solver": "Truncated CG",
            "Subproblem Model": "Kelley Sachs"
        }
    },
    "Status Test" :
    {
        "Constraint Tolerance":  1e-9,
        "Gradient Tolerance": 1e-7,
        "Iteration Limit": 200
    }
}
cs_dict = {**params_dict, **{"Step": {"Type": "Composite Step"}}}
ip_dict = {**params_dict, **{"Step": {"Type": "Interior Point"}}}
al_dict = {**params_dict, **{"Step": {"Type": "Augmented Lagrangian"}}}
fl_dict = {**params_dict, **{"Step": {"Type": "Fletcher"}}}
ls_dict = {**params_dict, **{"Step": {"Type": "Line Search"}}}
tr_dict = {**params_dict, **{"Step": {"Type": "Trust Region"}}}
my_dict = {**params_dict, **{"Step": {"Type": "Moreau-Yosida Penalty"}}}


def test_HS4_ls():
    print("HS4")
    run_HS4(ROL.StdVector, ls_dict)
    run_HS4(NumpyVector, ls_dict)

def test_HS4_tr():
    print("HS4")
    run_HS4(ROL.StdVector, params_dict)
    run_HS4(NumpyVector, params_dict)

def test_HS13_my():
    print("HS13")
    run_HS13(NumpyVector, my_dict)
    run_HS13(NumpyVector, my_dict)

# def test_HS28_al():
#     run_HS28(ROL.StdVector, params_dict)
#     run_HS28(NumpyVector, params_dict)

def test_HS28_cs():
    print("HS28")
    run_HS28(ROL.StdVector, cs_dict)
    run_HS28(NumpyVector, cs_dict)

# def test_HS28_my():
#     run_HS28(ROL.StdVector, my_dict)
#     run_HS28(NumpyVector, my_dict)

def test_HS29_ip():
    print("HS29")
    run_HS29(ROL.StdVector, ip_dict)
    run_HS29(NumpyVector, ip_dict)

def test_HS29_al():
    print("HS29")
    run_HS29(ROL.StdVector, al_dict)
    run_HS29(NumpyVector, al_dict)

def test_HS29_my():
    print("HS29")
    run_HS29(ROL.StdVector, my_dict)
    run_HS29(NumpyVector, my_dict)

# def test_HS29_fl():
#     run_HS29(ROL.StdVector, fl_dict)
#     run_HS29(NumpyVector, fl_dict)
