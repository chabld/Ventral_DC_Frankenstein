#Code to extract and merge: Brainstem Navigator's geniculate nuclei with minimal probability threshold, same as was done for the FreeSurfer hippocampal/brainstem/thalamic segmentations, JHU DTI atlas's cerebral peduncles, and CIT168's nuclei (overlapping with ventral-DC)
#Mostly shamelessly Claude-generated (with several tweaks) 
#It is intended to be run from a root directory that includes 3 subdirectories, i.e. "2b.diencephalicNucleiAtlas_MNI/" (including LG and MG volumes), "JHU/" (including the JHU atlas), and "CIT168/" including the CIT168toMNI152-2009c_det atlas (split using CIT168_splitter.py).

import nibabel as nib
import numpy as np
from scipy.ndimage import affine_transform

######################################## 

# --- MGN/LGN: Brainstem Navigator files ---
files = {
    1: f"2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/LG_l.nii.gz",
    2: f"2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/LG_r.nii.gz",
    3: f"2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/MG_l.nii.gz",
    4: f"2b.DiencephalicNucleiAtlas_MNI/labels_probabilistic/MG_r.nii.gz",
}

# --- load reference for affine/header ---
ref_img = nib.load(files[1])
# --- load all probability maps ---
probs = {label: nib.load(path).get_fdata() for label, path in files.items()}

# --- stack: index 0 = "background/none", 1-4 = LG_l, LG_r, MG_l, MG_r ---
label_ids = sorted(probs.keys())
stack = np.stack(
    [np.zeros_like(probs[label_ids[0]])] + [probs[i] for i in label_ids],
    axis=-1)

# --- mask out non-zeros ---
total = sum(probs.values())
valid = total >= 0.0001
# --- whichever label is highest probability in voxel is assigned ---
argmax_idx = np.argmax(stack, axis=-1)  # 0=none, 1=LG_l, 2=LG_r, 3=MG_l, 4=MG_r
labels = np.where(valid, argmax_idx, 0).astype(np.int16)

######################################## 

# --- Cerebral peduncles: JHU DTI atlas file ---
jhu_img = nib.load("JHU/JHU-ICBM-labels-1mm.nii.gz")  
jhu_data = jhu_img.get_fdata()

# Only interested in cerebral peduncle R and L
peduncle_r = (jhu_data == 15)
peduncle_l = (jhu_data == 16)
#assign new label IDs (5, 6) that don't collide with 1-4 already used
combined = labels.copy()
combined[peduncle_r] = 5  # cerebral peduncle R
combined[peduncle_l] = 6  # cerebral peduncle L
# note: if a voxel overlap with LGN/MGN, peduncle silently overwrites.
overlap = ((peduncle_r | peduncle_l) & (labels > 0)).sum()

######################################## 

# --- CIT168 atlas file ---
# /!\ this volume was splitted to separate left and right with the script "CTI168_splitter.py"
cit168_img = nib.load("CIT168/CIT168toMNI152-2009c_det_split.nii.gz")

# because it's MNI 152 2009c specifically, resample CIT168 onto the JHU/geniculate grid
# --- inputs ---
targ_img = jhu_img
src_data = cit168_img.get_fdata()
src_affine = cit168_img.affine
targ_affine = targ_img.affine
targ_shape = targ_img.shape

# --- compute voxel-to-voxel transform ---
# affine_transform maps OUTPUT voxel coords to INPUT voxel coords, so we need
# the inverse: src_affine^-1 @ targ_affine
src_affine_inv = np.linalg.inv(src_affine)
voxel_to_voxel = src_affine_inv @ targ_affine
matrix = voxel_to_voxel[:3, :3]
offset = voxel_to_voxel[:3, 3]

# --- resample with nearest-neighbor (order=0) - mandatory for integer labels ---
resampled = affine_transform(
    src_data,
    matrix=matrix,
    offset=offset,
    output_shape=targ_shape,
    order=0,          # nearest-neighbor, no interpolation of label values
    mode='constant',
    cval=0)

resampled = resampled.astype(np.int16)
out_img = nib.Nifti1Image(resampled, targ_affine, targ_img.header)
out_img.header.set_data_dtype(np.int16)

# structures to keep in the original CIT168 labels:
base_ids = {
    7: "Substantia_Nigra_pars_compacta",
    8: "Red_Nucleus",
    9: "Substantia_Nigra_pars_reticulata",
    10: "Parabrachial_Pigmented",
    11: "Ventral_Tegmental_Area",
    14: "Hypothalamus",
    16: "Subthalamic_Nucleus",
}

#in split volumes, left are original IDs and right ID + 16 
keep_values = list(base_ids.keys()) + [b + 16 for b in base_ids.keys()]
cit_filtered = np.where(np.isin(resampled, keep_values), resampled, 0)

# --- final label mapping ---
# left:  orig id (7-16)   -> +40  => 47-56
# right: orig id+16 (23-32) -> +44 => 67-76  (keeps same final scheme as before)
lookup = np.zeros(33, dtype=np.int16)  # covers label values 0-32
for b in base_ids:
    lookup[b] = b + 40          # left
    lookup[b + 16] = b + 60     # right (equivalent to +44 on the +16 id)

cit_shifted = lookup[cit_filtered.astype(int)]

# --- check overlap ---
overlap_final = ((cit_shifted > 0) & (combined > 0)).sum()
if overlap_final > 0:
    print(f"WARNING: {overlap_final} voxels overlap between CIT168 and geniculate/peduncle labels")
    print("Overlapping CIT168 (shifted) values:",
          np.unique(cit_shifted[(cit_shifted > 0) & (combined > 0)]))
    print("Overlapping combined (geniculate/peduncle) values:",
          np.unique(combined[(cit_shifted > 0) & (combined > 0)]))

# --- priority rules ---
# geniculate (1-4) wins over CIT168, which wins over JHU cerebral peduncle 
final = cit_shifted.copy()
genic_mask = np.isin(combined, [1, 2, 3, 4])
final[genic_mask] = combined[genic_mask]
peduncle_mask = np.isin(combined, [5, 6]) & (final == 0)
final[peduncle_mask] = combined[peduncle_mask]

# --- save ---
final_img = nib.Nifti1Image(final, targ_affine, targ_img.header)
final_img.header.set_data_dtype(np.int16)
nib.save(final_img, "VentralDC.MNI152.nii.gz")

