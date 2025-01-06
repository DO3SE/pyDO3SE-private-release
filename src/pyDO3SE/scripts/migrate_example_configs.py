from json.decoder import JSONDecodeError
import json
from pyDO3SE.version import config_version
from pyDO3SE.Config.config_migration import MigrationError, Migrations
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path + "/../")
print(sys.path)
EXAMPLE_DIRS = ["./examples", "./tests"]


if __name__ == "__main__":
    print("Migrating configs")
    # It also selects the first csv file in the directory that is not the external data file
    for EXAMPLE_DIR in EXAMPLE_DIRS:
        config_paths = [[
            os.path.basename(dr),  # Dir name
            dr,  # Directory path
            [fl for fl in f
             if fl[-5:] == '.json' and fl not in [
                 "config.json",
                 "final_state.json",
                 "processed_config.json",
                 "initial_state.json",
                 "expected_state.json",
                 "preprocess_map.json",
                 "variable_map.json",
                 "base_state",
             ]],  # Target file in dir
        ]
            for dr, d, f in os.walk(EXAMPLE_DIR)
            if os.path.isdir(dr)]
        for name, dir, files in config_paths:
            if dir == "outputs":
                continue
            if "/runs/" in dir:
                continue
            if "/e_state_overrides_field_maps/" in dir:
                continue
            for file in files:
                config_location = dir + "/" + file
                with open(config_location) as config_file_data:
                    try:
                        config = json.load(config_file_data)
                        input_version = config.get('VERSION', config_version + 1)
                        input_comment = config.get("COMMENT", "")
                        if not input_comment and not config.get('VERSION', 0):
                            print(f"Skipping config: {config_location} - No comment")
                            continue

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
                    except Exception as e:
                        print(f"Failed to migrate config: {config_location} - {e}")
                        if "-v" in sys.argv:
                            raise e
