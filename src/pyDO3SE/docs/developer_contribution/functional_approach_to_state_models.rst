###########################################################################
Approach to state based models(WIP) using a Functional Programming approach
###########################################################################
OUT OF DATE!

.. contents:: Contents


This is a work in progress documentation of good practice for writing
readable code that avoids common pitfalls such as rabbit hole variables.
The general aim is to be explicit on the inputs, internal state changes
and outputs of each model and sub model. It uses Functional Programming
approaches.


************
Key concepts
************

Config, Internal/Model State and External State
===============================================

The model data is broken down into config, state and data. - **Config**
is the configuration of the model such as number of runs and type of
process. It should be set at the start and be unchanged throughout the
model. - **State - (Internal)Model State** is the state inside the
running model such as time of day. Model State is modified by processes
within the model - **Data - External State** is inputed at the start of
the model such as recorded weather data. It should be set at the start
and be unchanged throughout the model

Functional Programming
======================

Immutability
------------

Functional Programming relies on 2 main concepts; **immutability** and
**Function Purity** to ensure Immutability Immutable variables
cannot be modified without creating a copy. Ensuring that variables are
immutable ensures that the value of that variable remains constant
throughout the function. This is important for readability of the code
and also for avoiding 'side effects' in your functions. It also allows
safe multithreading which can speed up simulations.

Example of Variables being modified:

.. code:: python

    a = 1
    print(a) # 1
    a = a + 1
    print(a) # 2

Functional Approach

.. code:: python

    a = 1
    print(a) # 1
    b = a + 1
    print(a) # 1
    print(b) # 2

In the second approach we can guarantee that anytime we call 'a' we get
the same value.

Pure Functions
~~~~~~~~~~~~~~

Pure functions are those that have no side effects. This means that
running the function will not have any effect on anything outside the
function and will not modify any of the inputs to the function. This
means that the function is deterministic, i.e given a set of inputs it
will always provide the same output. Ths makes testing and modulerising
of code much easier and safer.

Example of an impure function:

.. code:: python

    a = 1

    def double_func():
        global a # accessing external variable
        a = a * 2

    print(a) # 1
    double_func()
    print(a) # 2

Pure Version

.. code:: python

    a = 1

    def double_func(a_in):
        a_out = a_in * 2
        return a_out

    print(a) # 1
    new_a = double_func(a)
    print(a) # 1
    print(new_a) # 2

Use of Named Tuple Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the named tuple class allows us to predefine all of the input,
state and output varables into one place and makes them immutable. It
also allows for easier identification of where the state is being
modified.

**SEE model\_data.md**

Plugin Interfaces
-----------------

Abstracting functions enables better code modulerising. We can use
interfaces to act as a middleware between 2 different systems. In
practice this means that we can create an interface for any plugin we
integrate so that the main software can work with the plugin even if the
plugin uses a different data strategy.

**SEE plugin\_system.md**

Process Runners to update state
-------------------------------

Pure functions should not modify the state. Instead they should return a
copy of the state.

Each function should take the input state and output a modified copy of
that state without effecting the original state. The global state can be
then set to equal this new state. This ensures that we can easily test
each process and be explicit on what it is modifying.

**SEE process\_runner.md**

Testing
-------

Testing allows us to be certain that the outputs of the model are
predictable when we change code. I.e we need to know we havent broken an
unrelated part of the model when we modify a piece of code.

Snapshot Testing
~~~~~~~~~~~~~~~~

The large quantity of variables used in these models can make testing
messy. To quickly take a snapshot of the output variables of a method we
can use snapshot testing. Future test runs will then check the output
against this snapshot and fail if there are any changes.

Function Mocks in Tests
~~~~~~~~~~~~~~~~~~~~~~~

Often it is desirable to test a function without testing the
functionality of the embedded processes. For example testing the main
function calls the initialization functions without actually running the
initialization functions. We will be using unittest.mock which is
included with python >3.3
https://docs.python.org/dev/library/unittest.mock.html

Other useful coding guidelines
------------------------------

Linting
-------

Code should use a standard python linter to ensure consistency (E.g.
flake or pylint) This provides warnings when the code is 'messy'.

Type hints
----------

Type hints help us know the expected variable type:

.. code:: python

    variable_a = "foo" # no type hint
    variable_b: str = "bar" # with type hint

Model(Primary) Module Code Structure
====================================

For modules that are primary(i.e. cannot be swapped out for a different
model) Model module code is structured as follows:

Imports
-------

imports allow you to bring in constants and pure functions from external
modules. Where possible the inports should be named (Avoid use of \*)

.. code:: python


    from pyDO3SE.constants.physical_constants import T0, R, DRATIO
    from pyDO3SE.met.helpers import saturated_vapour_pressure
    from pyDO3SE.photosynthesis.constants import *
    from pyDO3SE.photosynthesis.helpers import calc_fO3_h, calc_f_LS

Input\_Shape
------------

This defines all the inputs to this part of the model. These should be
not changed throughout the model.

.. code:: python

    class Inputs(NamedTuple):
        """Inputs"""
        input_variable_a: int = None # Some info on the variable

State\_Shape
------------

This defines all the state variables within the model. These can be
altered over the lifetime of the model but should be copied for each
state change (see model run below). State is only required if the model
runs a loop within itself that will effect state on each loop.

.. code:: python

    class State(NamedTuple):
        """Variables that change over the model run"""
        state_variable_a: float = None # Some info about the variable

Output\_Shape
-------------

The output shape defines the variables that were outputed from this
module. These should only be variables generated from within the module.
As the input is immutable

.. code:: python

    class Output(NamedTuple):
        """Outputs"""
        output_variable_a: str = None # Some info about the variable

Methods
-------

Define methods used here. They should take the model input and current
state as inputs and output a copy of the state

.. code:: python

    def modify_state(inputs: Inputs, state_1: State) -> State:
        """ Process to modify state...

        State Effects:
        - state_var_01
        """
        # Define variables taken from inputs
        input_variable_a: int = inputs.input_variable_a
        state_variable_a: float = inputs.state_variable_a

        # Do some stuff to the state variable

        new_state_variable_a = input_variable_a * state_variable_a

        # Return a copy of the state
        return state_1._replace(
            state_variable_a=new_state_variable_a,
        )

Model Run Method
----------------

At the bottom of the file we should run the model using the data and
methods above then output using the Output Class. Note we can input
either by passing an instance of the input class or as arguments..

.. code:: python

    def run_model(inputs_in: Inputs = None, state_in: State = None, **kwargs) -> Output:
        """ Run model and return output tuple"""

        # gets input and state from method input
        inputs = inputs_in or Inputs(**kwargs)
        state_1 = state_in or State()

        # get modified state variables
        state_2 = state_1._replace(state_variable_a=1.1)
        state_3 = modify_state(inputs, state_2)

        output = Output(
            output_variable_a=state_3["state_variable_a"],
        )

        return Output()

Running the model
-----------------

.. code:: python

    from module import run_model, Inputs

    #1. Using the input class
    inputs = Inputs(
        input_variable_a=1
    )
    out_a = run_model(inputs_in=inputs)
    assert isinstance(out_a, Output)
    print(tuple(out_a))
    # >>> Output(output_variable_a = 2.2)

    #2. Using named arguments
    out_ab = run_model(input_variable_a=1)
    assert isinstance(out_b, Output)
    print(tuple(out_b))
    # >>> Output(output_variable_a = 2.2)

Model(Plugin) Module Code Structure
===================================

For models that can be classed as plugins (i.e they could be swapped out
for a different model).

These will be the same as primary modules except all the plugins will
share the same input and output

Example
-------

.. code:: python

    class Model_Input_shape(NamedTuple):
        foo: int
        bar: int

    class Model_Output_Shape(NamedTuple):
        result: int

    def plugin_model_a(input: Model_Input_shape) -> Model_Output_Shape:
        foo = input.foo
        bar = input.bar
        result = foo * bar
        return Model_Output_Shape(
            result=result,
        )

    def plugin_model_b(input: Model_Input_shape) -> Model_Output_Shape:
        foo = input.foo
        bar = input.bar
        result = foo + bar
        return Model_Output_Shape(
            result=result,
        )

    input = Model_Input_shape(
        foo = 2,
        bar = 3,
    )

    plugin_model_a(input)
    # returns (result: 6) i.e 2 * 3

    plugin_model_b(input):
    # returns (result: 5) i.e 2 + 3
