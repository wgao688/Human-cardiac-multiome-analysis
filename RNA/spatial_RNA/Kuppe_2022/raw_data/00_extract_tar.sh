#!/bin/bash

for tarball in $PWD/*.tar.gz; do
	tar -xvf $tarball
done
