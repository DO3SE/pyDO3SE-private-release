import sys
from json.decoder import JSONDecodeError
import os
import json
from pyDO3SE.version import config_version
from pyDO3SE.Config.config_migration import MigrationError, Migrations
EXAMPLE_DIR = "./examples"


if __name__ == "__main__":
    print("Migrating configs")
    # It also selects the first csv file in the directory that is not the external data file
    config_paths = [[
        os.path.basename(dr),  # Dir name
        dr,  # Directory path
        [fl for fl in f
         if fl[-5:] == '.json' and fl not in ["config.json", "final_state.json", "processed_config.json"]],  # Target file in dir
    ]
        for dr, d, f in os.walk(EXAMPLE_DIR)
        if os.path.isdir(dr)]
    for name, dir, files in config_paths:
        for file in files:
            config_location = dir + "/" + file
            with open(config_location) as config_file_data:
                try:
                    config = json.load(config_file_data)
                    input_version = config.get('VERSION', 0)
                    migrated_config = Migrations.run_migrations(
                        config, input_version, verbose=False)

                    out_dir = os.path.dirname(config_location)
                    out_file_name = os.path.splitext(os.path.basename(config_location))[
                        0] + f"_{config_version}.json"
                    with open(config_location, 'w') as out_location:
                        json.dump(migrated_config, out_location, indent=4)

                except MigrationError as e:
                    print(f"Failed to migrate config: {config_location} - {e.message}")
                    if "-v" in sys.argv:
                        raise e
                except JSONDecodeError as e:
                    print(f"Failed to load file: {config_location}")
                    if "-v" in sys.argv:
                        raise e
