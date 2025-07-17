from pyDO3SE.settings import settings
from pyDO3SE.Output.Output_Shape import output_fields

def get_multi_dimension_output_fields(fieldId: str):
    field_data = next(o for o in output_fields if o.id == fieldId)
    field_additional_multi_dimension_field_ids = []
    if field_data.per_iL:
        field_additional_multi_dimension_field_ids += [
            f"{fieldId}_iL_{iL}" for iL in range(settings().MAX_NUM_OF_CANOPY_LAYERS)
        ]
    else:
        field_additional_multi_dimension_field_ids += [fieldId]
    if field_data.per_iP:
        field_additional_multi_dimension_field_ids = [
            f"{ff}_iP_{iP}"
            for ff in field_additional_multi_dimension_field_ids
            for iP in range(settings().MAX_NUM_OF_LEAF_POPULATIONS)
        ]
    if field_data.per_iCH:
        field_additional_multi_dimension_field_ids += [
            f"{fieldId}_iCH_{iCH}" for iCH in range(settings().MAX_NUM_OF_CUSTOM_LAYERS)
        ]
    return field_additional_multi_dimension_field_ids
