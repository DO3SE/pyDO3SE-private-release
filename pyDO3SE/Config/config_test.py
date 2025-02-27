# import pytest

# from pyDO3SE.Config import Config_Shape, init_config, get_config, deinit_config
# from pyDO3SE.Config.Config import __config  # noqa: F401


# DEMO_CONFIG = Config_Shape()


# def test_init_config():
#     ''' Check we can initialise the config'''
#     with pytest.raises(Exception):
#         get_config()
#     init_config(DEMO_CONFIG)
#     assigned_config = get_config()
#     assert assigned_config == DEMO_CONFIG

#     # deinitializes config
#     deinit_config()
#     with pytest.raises(Exception):
#         get_config()

#     # raises exception if trying to initialize config twice
#     with pytest.raises(Exception):
#         init_config(DEMO_CONFIG)
#         init_config(DEMO_CONFIG)
#     deinit_config()

#     # raises exception if trying to deinitialize config before initialized
#     with pytest.raises(Exception):
#         deinit_config()


# def test_config_is_private():
#     ''' check that modifiying the __config object does not change the get_config output'''
#     global __config
#     assert __config is None
#     __config = DEMO_CONFIG

#     with pytest.raises(Exception):
#         get_config()  # noqa: F841


# def test_cannot_pass_invalid_config_object_to_init():
#     with pytest.raises(Exception):
#         get_config()
#     invalid_config = {'a': 'foo'}
#     with pytest.raises(Exception):
#         init_config(invalid_config)

#     invalid_config_b = Config_Shape(invalid_config)
#     with pytest.raises(Exception):
#         init_config(invalid_config_b)

#     with pytest.raises(Exception):
#         get_config()
