from .Output_Shape import default_output_fields, output_fields_map


class TestOutputFieldMap:

    def test_should_have_all_default_output_fields_in_output_fields(self):
        for f in default_output_fields:
            try:
                field = output_fields_map[f]
            except KeyError:
                raise KeyError(
                    f"{f} not in output fields map. Run cli command `pyDO3SE_cli available-outputs")
