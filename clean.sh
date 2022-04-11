#!/bin/bash

remove_dir() {
    find . -type d -name $1 \
        -exec echo "removing {}" \; \
        -exec rm -dr {} +
}

remove_ext() {
    find . -type f -name "*$1" \
        -exec echo "removing {}" \; \
        -exec rm {} +
}

remove_dir ".ipynb_checkpoints"
remove_dir "__pycache__"
remove_dir ".vscode"

for EXT in '.aux' '.bbl' '.blg' '.fdb_latexmk' '.fls' '.log' '.out' '.synctex.gz' '.toc'; do
    remove_ext $EXT
done
