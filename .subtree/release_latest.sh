#!/bin/bash
set -e

.subtree/clone_and_pull.sh pull RELEASE
.subtree/snapshot_versions.sh
