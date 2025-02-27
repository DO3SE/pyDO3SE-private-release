# %%
# TESTING PLF_Value
from math import isclose
from typing import List


def PLF_value(points: List[List[float]], x: float) -> float:
    n = len(points[1]) - 1
    print(n)
    if x < points[0][0]:
        print("a")
        y = points[1][0]
    elif x >= points[0][n]:
        print("b")
        y = points[1][n]
    else:
        print("c")
        bx = points[0][0]
        by = points[1][0]
        print(bx, by)

        for i in range(1, n + 1):
            ax = bx
            ay = by
            bx = points[0][i]
            by = points[1][i]

            # Skip zero-width pieces(this should be equivalent to an # equality check, but checking floating point equality is evil # and the compiler warns about it)
            if isclose(abs(ax - bx), 0.0):
                continue

            if (x <= bx):
                y = ay + (by - ay) * ((x - ax) / (bx - ax))
                break
    return y


def run_plf_value(dd): return PLF_value([[0.0, 0.0, 0.0, 0.0, 0.0, 47.0, 47.0, 47.0, 92.0, 92.0], [
    0.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.0]], dd)


# assert run_plf_value(dd=364) == 0.0 # ddadj = 364
assert run_plf_value(dd=118 - 118) == 1.0  # ddadj = 0 (365)
# assert run_plf_value(dd=165 - 118) == 1.0
# assert run_plf_value(dd=166 - 118) == 0.98
# TODO: for dd = 210 we are getting 0.0 instead of 0.1 (Value from UI)
assert run_plf_value(dd=210 - 118) == 0.1
# assert run_plf_value(dd=211 - 118) == 0.0
