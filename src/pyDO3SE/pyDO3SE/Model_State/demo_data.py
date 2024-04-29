from .Model_State import Model_State_Shape
example_state = Model_State_Shape()
example_state.temporal.td = 0
example_state.canopy_component[0].td = 0
example_state.prev_hour = Model_State_Shape()
example_state.prev_hour.temporal.td = 0
example_state.prev_hour.canopy_component[0].td = 0
