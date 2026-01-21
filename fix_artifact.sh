#!/bin/bash

# Usage: ./fix_artifact.sh <HASH> <URL>

HASH=$1
URL=$2
ARTIFACT_DIR="$HOME/.julia/artifacts/$HASH"

if [ -z "$HASH" ] || [ -z "$URL" ]; then
    echo "Usage: ./fix_artifact.sh <HASH> <URL>"
    echo "Example: ./fix_artifact.sh e834... https://github.com/..."
    exit 1
fi

echo "------------------------------------------------"
echo "Fixing Julia Artifact..."
echo "Hash: $HASH"
echo "URL:  $URL"
echo "------------------------------------------------"

# 1. Create the directory
if [ -d "$ARTIFACT_DIR" ]; then
    echo "[!] Directory exists. Removing partial/corrupt files..."
    rm -rf "$ARTIFACT_DIR"
fi
mkdir -p "$ARTIFACT_DIR"

# 2. Download
echo "[*] Downloading tarball..."
wget -q --show-progress -O /tmp/temp_artifact.tar.gz "$URL"

if [ $? -ne 0 ]; then
    echo "[X] Download failed! Check your URL or internet connection."
    exit 1
fi

# 3. Extract
echo "[*] Extracting..."
tar -xzf /tmp/temp_artifact.tar.gz -C "$ARTIFACT_DIR"

if [ $? -eq 0 ]; then
    echo "[✓] Success! Artifact installed to $ARTIFACT_DIR"
    echo "    You can now return to Julia and retry Pkg.add()"
    rm /tmp/temp_artifact.tar.gz
else
    echo "[X] Extraction failed."
fi
