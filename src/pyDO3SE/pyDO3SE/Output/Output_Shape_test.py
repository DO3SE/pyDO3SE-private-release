from .Output_Shape import default_output_fields, output_fields_map


class TestOutputFieldMap:

    def test_should_have_all_default_output_fields_in_output_fields(self):
        for f in default_output_fields:
            if f in ["lai", "lai_brown", "layer_lai_brown"]:
                # Skip because these are not in output_fields_map
                continue
            try:
                field = output_fields_map[f]
            except KeyError:
                raise KeyError(
                    f"{f} not in output fields map. Run cli command `pyDO3SE config available-outputs`")
