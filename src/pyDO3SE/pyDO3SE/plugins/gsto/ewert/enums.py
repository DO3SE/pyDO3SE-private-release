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



class AdjustNegativeAnMethods(Enum):
    """Methods for adjusting negative An values.


    FALSE = "false"
        Do not adjust negative An values
    ALLOW = "allow"
        ALlow negative An values
    LAST_RESORT = "last_resort"
        Allow negative An values as a last resort in cubic solver
    CLIP = "clip"
        Clip negative An values to 0
    CLIP_IN_LOOP = "clip_in_loop"
        Clip negative An values to 0 in the An loop

    """

    FALSE = "false"
    ALLOW = "allow"
    LAST_RESORT = "last_resort"
    CLIP = "clip"
    CLIP_IN_LOOP = "clip_in_loop"