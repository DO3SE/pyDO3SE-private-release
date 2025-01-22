# DO3SE Open Repository

This is a repository for hosting released versions of the DO3SE model.

It uses a Git SubTree workflow to pull the release version of the sub repositories.

DO NOT MODIFY THE CONTENTS OF THE 'src' DIRECTORY DIRECTLY.

## Updating versions

To pull changes to the sub repositories, use the following command:

```bash
./.subtree/clone_and_pull.sh pull main
```

This will pull the latest changes from the 'main' branch of the sub repositories.

To snapshot a version run the following command:

```bash
./.subtree/snapshot_versions.sh
```

This will create a new commit with the current versions of the sub repositories.

If a sub repository has changed and the version of pyDO3SE hasn't also been updated
the command will throw an error. This is to ensure that we can keep the version of
pyDO3SE-open in sync with the version of pyDO3SE.

To pull specific versions of the sub repositories, use the following command:

```bash
./.subtree/pull_specific_versions.sh <pyDO3SE_VERSION> <THERMAL_TIME_VERSION> <do3se_phenology_VERSION> <do3se_met_VERSION>
```

## Environment setup

You will need to first install `uv` see the [uv documentation](https://docs.astral.sh/uv/getting-started/installation/) for more information.

If setting up for the first time, run the following commands:

- `scripts/env.init.sh`

If updating the environment, run the following commands:

- `scripts/env.refresh.sh`
