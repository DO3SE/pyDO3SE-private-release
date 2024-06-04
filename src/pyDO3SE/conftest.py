def pytest_addoption(parser):
    parser.addoption("--arg", action="store", default="")
    parser.addoption("--skip-outputs", action="store", default="False")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.arg
    if 'arg' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("arg", [option_value])
