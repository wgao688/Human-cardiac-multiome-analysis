#!/bin/bash

tarball_directory="tarballs"
mkdir -p $tarball_directory

# for each tarball in the current directory, extract it and then move the tarball into the tarball directory
for tarball in *.tar.gz; do
	tar_name=$(basename "$tarball" .tar.gz)

	# extract the tarball
	tar -xvf  "$tarball"
	# rename the GeneFull_Ex50pAS/ dir to the tarball name
	mv "GeneFull_Ex50pAS/" $tar_name

	# move the tarball into the tarball directory
	mv $tarball $tarball_directory
done
