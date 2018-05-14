import copy
from hs.hs4 import run_HS4
from hs.hs13 import run_HS13
from hs.hs28 import run_HS28
from hs.hs29 import run_HS29
import ROL
from ROL.numpy_vector import NumpyVector
params_dict = {
    "Step": {
        "Type": "Moreau-Yosida Penalty",
        "Moreau-Yosida Penalty":  {
            "Initial Penalty Parameter": 1.,
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

def get_step_dict(steptype):
    dict_ = copy.deepcopy(params_dict)
    dict_["Step"]["Type"] = steptype
    return dict_

cs_dict = get_step_dict("Composite Step")
ip_dict = get_step_dict("Interior Point")
al_dict = get_step_dict("Augmented Lagrangian")
fl_dict = get_step_dict("Fletcher")
ls_dict = get_step_dict("Line Search")
tr_dict = get_step_dict("Trust Region")
my_dict = get_step_dict("Moreau-Yosida Penalty")


def test_HS4_ls():
    print("HS4")
    run_HS4(ROL.StdVector, ls_dict)
    run_HS4(NumpyVector, ls_dict)

def test_HS4_tr():
    print("HS4")
    run_HS4(ROL.StdVector, tr_dict)
    run_HS4(NumpyVector, tr_dict)

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
