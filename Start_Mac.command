#!/bin/bash
# CLAVAGEN - double-click to start (macOS). The CLAVAGEN.app icon does the same without a Terminal window.
cd "$(dirname "$0")" || exit 1
exec ./tools/launch_mac.sh
