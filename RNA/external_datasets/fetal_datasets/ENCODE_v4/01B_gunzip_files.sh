#!/bin/bash

for file in ENC*/raw*/*; do
	gzip $file
done
