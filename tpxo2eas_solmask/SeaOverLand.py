from __future__ import print_function, division
import sys
import numpy as np
import numpy.ma as ma
import netCDF4 as NC
###############################
# INPUT
# Tpxo file
infile_pathname=sys.argv[1]
infile_var='tide_z'
# In/Out files Dimensions names 
lat_idx='y' #'nav_lat'
lon_idx='x' #'nav_lon'
time_idx='time_counter'
# Tpxo Mask file
mask_pathname=sys.argv[2]
mask_var='mask4tpxo'
# NEMO Mask file
nemo_mask=sys.argv[3]
nemo_tmask='tmask'
# loop num
sol_loop=int(sys.argv[5])
#
# OUTPUT
# Output file (MUST contain the same dimensions as the input one)
outfile_pathname=sys.argv[4]
outfile_var='tide_z'
###############################

# Open input tpxo files to be modified 
infile = NC.Dataset(infile_pathname,'r')
infield=np.array(infile.variables[infile_var][:])

# Open the mask 
mask = NC.Dataset(mask_pathname,'r')
inmask=np.array(mask.variables[mask_var][:])

# Inizialize the new field to the old one
#new_field=np.array(infield)

# Close the inputs nc 
infile.close()
mask.close()

# Open the file to write the new field
outfield = NC.Dataset(outfile_pathname,'r+')
# Build the new field
tide_z_masked=outfield.createVariable(outfile_var, 'f4', (time_idx, lat_idx , lon_idx,))
tide_z_masked.units = 'm'

# Mask the field
#mult_field = inmask[:]*infield[:]
#input_matrix=mult_field

threshold = 0
mask2interp = np.abs(inmask) == threshold
masked_field = np.ma.masked_where(mask2interp,infield)
input_matrix=masked_field

# Run the SeaOverLand
nloop=sol_loop

infill_value = input_matrix.fill_value
#print ('infill_value ',infill_value)

output_matrix = ma.copy(input_matrix)  # using ma.copy to prevent future field modifications
if np.sum(output_matrix.mask) == 0:  # nothing to fill
        print('WARNING. Field does not have any land points or wrong field type. Exiting.', file=sys.stderr)
else:
        # iterations loop
        print(nloop)
        for loop in range(nloop):
            # Create a nD x 8 matrix in which, last dimension fixed, the other dimensions
            #  contains values that are shifted in one of the 8 possible directions
            # of the last two axes compared to the original matrix
            shift_matrix = ma.array(np.empty(shape=((8,) + output_matrix.shape)),
                                    mask=True, fill_value=infill_value, dtype=float)
            # up shift
            shift_matrix[0, ..., : -1, :] = output_matrix[..., 1:, :]
            # down shift
            shift_matrix[1, ..., 1:, :] = output_matrix[..., : -1, :]
            # left shift
            shift_matrix[2, ..., :, : -1] = output_matrix[..., :, 1:]
            # right shift
            shift_matrix[3, ..., :, 1:] = output_matrix[..., :, : -1]
            # up-left shift
            shift_matrix[4, ..., : -1, : -1] = output_matrix[..., 1:, 1:]
            # up-right shift
            shift_matrix[5, ..., : -1, 1:] = output_matrix[..., 1:, : -1]
            # down-left shift
            shift_matrix[6, ..., 1:, : -1] = output_matrix[..., : -1, 1:]
            # down-right shift
            shift_matrix[7, ..., 1:, 1:] = output_matrix[..., : -1, : -1]
            # Mediate the shift matrix among the third dimension
            mean_matrix = ma.mean(shift_matrix, 0)
            # Replace input masked points with new ones belonging to the mean matrix
            output_matrix = ma.array(np.where(mean_matrix.mask + output_matrix.mask, mean_matrix, output_matrix),
                                     mask=mean_matrix.mask, fill_value=infill_value, dtype=float)
            output_matrix = ma.masked_where(mean_matrix.mask, output_matrix)
            if np.sum(output_matrix.mask) == 0:  # nothing more to flood
                print('WARNING. Field does not have anymore land points,', str(loop + 1),
                      'steps were sufficient to flood it completely.', file=sys.stderr)
                break
# Open NEMO mesh_mask and apply the NEMO model tmask
meshmask = NC.Dataset(nemo_mask,'r')
nemo_land=meshmask.variables[nemo_tmask][0,0,:,:]

output_matrix=np.nan_to_num(output_matrix)
output_matrix=np.where(output_matrix==infill_value,0,output_matrix)

for time_output in range(0,len(output_matrix[:,0,0])):
    output_matrix[time_output,:,:]=output_matrix[time_output,:,:]*nemo_land[:,:]

# Write the new field to the out file 
tide_z_masked[:]=output_matrix[:] #*nemo_land[:]

# Close the outfile
outfield.close()
