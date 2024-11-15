from enum import Enum

class EwertLoopMethods(Enum):
    """We have two methods of converging on a A_n and gsto value.


    ITERATIVE = "iterative"
        Iteratively solve for A_n and gsto values
    CUBIC = "cubic"
        Use a cubic equation solver to solve for A_n and gsto values
        This is the faster method and is recommended for most cases.

    """

    ITERATIVE = "iterative"
    CUBIC = "cubic"