# SETUP INTERACTIVE WIDGETS
# @title Setup interactive
# @markdown code hidden

from functools import wraps
import ipywidgets as widgets
from IPython.display import display


def reset_button(defaults={}):
    def on_button_clicked(_):
        for k, v in defaults.items():
            k.value = v
    button = widgets.Button(description='Reset')
    button.on_click(on_button_clicked)
    display(button)


def add_label(w, label):
    return widgets.HBox([widgets.Label(label), w], layout={'width': '800px'})


def interact_plus_reset(**_widgets):
    default_vals = {wid: wid._trait_values['value'] for k, wid in _widgets.items()}
    reset_button(defaults=default_vals)

    def wrap(func):
        abox = widgets.VBox(list([add_label(v, k)
                            for k, v in _widgets.items()]), layout={'width': 'max-content'})

        @wraps(func)
        def inner(*args, **kwargs):
            return func(**args, **kwargs)
        display(abox, widgets.interactive_output(func, _widgets))
        return inner
    return wrap

# EXAMPLE
# @interact_plus_reset(
#   foo_with_really_really_long_name=widgets.FloatSlider(min=0, max=1, value=0.3, step=0.1),
#   bar=widgets.FloatSlider(min=0, max=1, value=0.3, step=0.1),
# )
# def get_leaf_f_phen(
#     foo_with_really_really_long_name,
#     bar,
#   ):
#   print(foo_with_really_really_long_name, bar)
