# %%
# Requires openpyxl

# %%
import pandas as pd
from pathlib import Path

column_map = {
    "DOY": "dd",
    "h": "hr",
    "AirT": "Ts_C",
    "RH": "RH",
    "Wind": "u",
    "[O3]": "O3",
    "Rain": "precip",
    "RadGlob": "Rn",
    "PPFD": "PPFD",
}

def correct_humidity(df):
    df['RH'] = df['RH'] / 100
    return df

def pull_data_from_agmip_template(
    file_path: Path,
    O3_Trt: str,
) -> pd.DataFrame:
    df_meteo = pd.read_excel(file_path, sheet_name='Meteorological data')
    df_gas = pd.read_excel(file_path, sheet_name = "Gas Concentration data")
    df_gas_filtered = df_gas[df_gas['O3_Trt'] == O3_Trt]
    df_gas_filtered = df_gas[['[O3]']]
    df_gas_filtered.reset_index()
    df = pd.merge(df_meteo, df_gas_filtered, left_index=True, right_index=True, suffixes=(False, False))
    df.columns = [column_map.get(c, c) for c in df.columns]
    df = df[list(column_map.values())]
    df = correct_humidity(df)

    return df

# Example
# data_file = "examples/xiaoji_2008/Xiaoji_2008.xlsx"
# df = pull_data_from_agmip_template(data_file, 'A-O3')
# df

import sys
if __name__=="__main__":
    args = sys.argv[1:]
    df = pull_data_from_agmip_template(*args[0:2])
    out = args[2]
    df.to_csv(out, index=False)
