# Use tar -xvf to unzip the .tar file. This produces several sub .tar.gz files. The script here creates a directory based on each of these .tar.gz files (without the .tar.gz extension), and then it moves the files into their respective directories. From there, it goes into each of these directories and uses tar -xvf to unzip the files, to generate the matrices

import os
import shutil
import subprocess

files_dir = os.getcwd() + "/"

# iterate over each file in the directory
for filename in os.listdir(files_dir):
    if filename.endswith(".tar.gz"):
        # remove the .tar.gz extension to create directory name
        directory_name = filename.replace(".tar.gz", "")

        directory_path = os.path.join(files_dir, directory_name)
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)

        # move the .tar.gz file into the directory
        src = os.path.join(files_dir, filename)
        dst = os.path.join(directory_path, filename)
        shutil.move(src, dst)

        # unzip the file using tar command
        subprocess.run(["tar", "-xvf", dst, "-C", directory_path])

print("Extraction completed.")
