#Claude-generated script to split the CIT168 atlas (using the X axis midline) and rearrange labels so left and right ROIs are distinguished 
#Used on "CIT168_Reinf_Learn_v1.1.0/MNI152-Nonlin-Asym-2009c/CIT168toMNI152-2009c_det.nii" (https://doi.org/10.17605/OSF.IO/R2HVK)
#This was visually checked with FreeSurfer and only the mammillary nucleus has an ambiguous split, with a one-voxel-thick column being potentially part of either left or right side

import nibabel as nib
import numpy as np

# --- input: original unresampled CIT168 det atlas, MNI152-2009c space ---
img = nib.load("CIT168toMNI152-2009c_det.nii.gz")
data = img.get_fdata()
affine = img.affine
shape = data.shape

# --- structure names and labels ---
names = {
    1: "Putamen",
    2: "Caudate",
    3: "Nucleus_Accumbens",
    4: "Extended_Amygdala",
    5: "Globus_Pallidus_externa",
    6: "Globus_Pallidus_interna",
    7: "Substantia_Nigra_pars_compacta",
    8: "Red_Nucleus",
    9: "Substantia_Nigra_pars_reticulata",
    10: "Parabrachial_Pigmented",
    11: "Ventral_Tegmental_Area",
    12: "Ventral_Pallidum",
    13: "Habenular_nucleus",
    14: "Hypothalamus",
    15: "Mammillary_Nucleus",
    16: "Subthalamic_Nucleus",
}

# --- compute world-space x-coordinate for every voxel ---
# MNI convention: x < 0 = left hemisphere, x >= 0 = right hemisphere
i, j, k = np.meshgrid(
    np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing='ij'
)
voxel_coords = np.stack([i, j, k, np.ones_like(i)], axis=-1)
world_coords = voxel_coords @ affine.T
world_x = world_coords[..., 0]

is_left = world_x < 0
is_right = world_x >= 0

# --- split: left keeps original label (1-16), right gets +16 (17-32) ---
split = np.zeros_like(data, dtype=np.int16)
left_mask = is_left & (data > 0)
right_mask = is_right & (data > 0)

split[left_mask] = data[left_mask].astype(np.int16)
split[right_mask] = data[right_mask].astype(np.int16) + 16

# --- save ---
out_img = nib.Nifti1Image(split, affine, img.header)
out_img.header.set_data_dtype(np.int16)
nib.save(out_img, "CIT168toMNI152-2009c_det_split.nii.gz")

# --- write matching labels.txt ---
with open("CIT168toMNI152-2009c_det_split_labels.txt", "w") as f:
    for orig_id, name in names.items():
        f.write(f"{orig_id},{name}_L\n")
    for orig_id, name in names.items():
        f.write(f"{orig_id + 16},{name}_R\n")
