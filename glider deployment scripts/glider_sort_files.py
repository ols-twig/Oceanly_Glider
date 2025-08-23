# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 09:29:52 2025

@author: OllieTwigge
"""

from glob import glob

mainpath = r'G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Solomon Islands 2025\Data\GLIDER_DATA\rawnc'
testpath = mainpath+r"\*.dbd.nc"

fullpath= list(glob(testpath, recursive=False))
#%%

import os
import shutil

mainpath = r'C:\temp\solo_25_glider_process\rawnc'

# Set your source and destination directories
src_dir = mainpath
dst_dir = r'C:\temp\solo_25_glider_process\rawnc_srt'
cmb_fil = r'C:\temp\solo_25_glider_process\rawnc_cmb'

# Make sure destination exists
os.makedirs(dst_dir, exist_ok=True)

# Loop through all files in source
for filename in os.listdir(src_dir):
    if not filename.endswith(".nc"):  # only process nc files
        continue

    # Get the first 4 digits
    prefix = filename[:4]

    # Create a subdirectory for that prefix if it doesn’t exist
    target_dir = os.path.join(dst_dir, prefix)
    os.makedirs(target_dir, exist_ok=True)

    # Move the file into the new directory
    src_path = os.path.join(src_dir, filename)
    dst_path = os.path.join(target_dir, filename)

    shutil.move(src_path, dst_path)
    print(f"Moved {filename} -> {target_dir}")

#%%
import os
import shutil

dst_dir = r'C:\temp\solo_25_glider_process\rawnc_srt'
cmb_fil = r'C:\temp\solo_25_glider_process\rawnc_cmb'

# The parent folder that contains the subfolders (0125, 0126, etc.)
src_parent = dst_dir
# The new folder where everything will be collected
dst_dir = cmb_fil
os.makedirs(dst_dir, exist_ok=True)

# Mapping of the original filenames to their desired suffix
file_map = {
    "unit_1104rawdbd.nc": ".dbd.nc",
    "unit_1104rawebd.nc": ".ebd.nc"
}

# Loop through each prefix folder
for prefix in os.listdir(src_parent):
    prefix_path = os.path.join(src_parent, prefix)
    if not os.path.isdir(prefix_path):
        continue

    # Loop through files in the subfolder
    for filename in os.listdir(prefix_path):
        if filename in file_map:
            src_path = os.path.join(prefix_path, filename)

            # Build new filename: prefix + mapped extension
            new_name = f"{prefix}{file_map[filename]}"
            dst_path = os.path.join(dst_dir, new_name)

            shutil.copy2(src_path, dst_path)
            print(f"Copied {src_path} -> {dst_path}")