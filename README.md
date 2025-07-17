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
