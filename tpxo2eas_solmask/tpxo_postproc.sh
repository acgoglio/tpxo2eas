#!/bin/bash
#set -u
set -e
#set -x 
###########################
# Year to process
YEAR2PP=2023
# Loop num in SeaOverLand procedure (at least 2 are required)
LOOPNUM=2

# Inputs
# old dataset to be proessed:
OLD_TPXOEXTR_PATH="/data/opa/mfs/Med_static/MFS_TPXO_V0/${YEAR2PP}/"
OLD_TPXOEXTR_PRE="tpxoextr" 
# mesh mask
MESHMASK="/work/oda/ag15419/PHYSW24_DATA/RIVERS/NEMO_DATA0_EAS6_PO/mesh_mask.nc"
# tpxo mask
TPXOMASK="/work/oda/ag15419/PHYSW24_DATA/DETIDING/mask4tpxo.nc"

# Output
# workdir and archive path
NEW_TPXOEXTR_PATH="/work/oda/ag15419/tmp/new_tpxoextr/"
NEW_TPXOEXTR_PRE=${OLD_TPXOEXTR_PRE}

###########################

# Workir and environment
if [[ -d $NEW_TPXOEXTR_PATH ]]; then
   OUTDIR="${NEW_TPXOEXTR_PATH}/${YEAR2PP}/"
   if [[ -d ${OUTDIR} ]]; then
      cp SeaOverLand.py ${OUTDIR}/
      cd $NEW_TPXOEXTR_PATH
      cd ${OUTDIR}/
   else
      mkdir ${OUTDIR}
      cp SeaOverLand.py ${OUTDIR}/
      cd $NEW_TPXOEXTR_PATH/
      cd ${OUTDIR}/
   fi
   echo "I am in Workdir: $(pwd)"
else
   echo "ERROR: Workdir $NEW_TPXOEXTR_PATH NOT found!"
   exit
fi

# cut the old field and create the new outfile templates 
module purge
module load anaconda
module load curl/7.70.0 cmake/3.17.3 gams/28.2.0 gcc_9.1.0/9.1.0 gcc_9.1.0/gempack/12.885 gcc_9.1.0/OpenBLAS/0.3.9 gcc_9.1.0/papi/6.0.0 gcc_9.1.0/R/3.6.1 modules mysql/5.7.28 ncl/6.6.2 sqlite/3.32.2 subversion/1.14.0 wgrib/1.8.1.0b impi20.1/19.7.217 impi20.1/esmf/8.0.1-intelmpi-64-g impi20.1/hdf5/1.12.0 impi20.1/hdf5-threadsafe/1.12.0 impi20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/papi/6.0.0 impi20.1/parallel-netcdf/1.12.1 impi20.1/petsc/3.13.2 impi20.1/zoltan/3.8 intel20.1/20.1.217 intel20.1/advisor intel20.1/boost/1.73.0 intel20.1/cdo/1.9.8 intel20.1/cnvgrib/3.1.1 intel20.1/eccodes/2.17.0 intel20.1/esmf/8.0.1-mpiuni-64-g intel20.1/esmf/8.0.1-mpiuni-64-O intel20.1/exactextract/545f0d6 intel20.1/g2lib/3.1.0 intel20.1/gdal/3.1.0 intel20.1/hdf5/1.12.0 intel20.1/hdf5-threadsafe/1.12.0 intel20.1/inspector intel20.1/itac intel20.1/libemos/4.5.9 intel20.1/libemos/4.5.9 intel20.1/magics/3.3.1 intel20.1/nco/4.9.3 intel20.1/ncview/2.1.8 intel20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/proj/7.0.1 intel20.1/R/4.0.2 intel20.1/szip/2.1.1 intel20.1/udunits/2.2.26 intel20.1/valgrind/3.16.0 intel20.1/vtune intel20.1/w3lib/2.0.6 intel20.1/wgrib2/2.0.8


# Loop on months
for MONTHLY_RUN in 01 02 03 04 05 06 07 08 09 10 11 12; do
   # Loop on days
   for TOWORKON in $( ls ${OLD_TPXOEXTR_PATH}/${OLD_TPXOEXTR_PRE}_${YEAR2PP}${MONTHLY_RUN}*.nc ); do 

      DATEON=$( echo $TOWORKON | cut -f 5 -d"_" | cut -f 1 -d".") 
      echo " I am working on date: $DATEON" 
      echo "Old file: $TOWORKON"

      # create the outfile
      TOOUT="${NEW_TPXOEXTR_PRE}_${DATEON}.nc"
      echo "New file ${TOOUT}"
      ln -sf $TOWORKON tmp_${DATEON}.nc
      ncks -x -v tide_z tmp_${DATEON}.nc ${TOOUT}
      rm tmp_${DATEON}.nc     
   done
done

echo "New outfile templates are ready!"

# Build the new field
module purge
module load anaconda
source activate mappyenv

# Loop on months
for MONTHLY_RUN in 01 02 03 04 05 06 07 08 09 10 11 12; do
   # Loop on days
   for TOWORKON in $( ls ${OLD_TPXOEXTR_PATH}/${OLD_TPXOEXTR_PRE}_${YEAR2PP}${MONTHLY_RUN}*.nc ); do

      DATEON=$( echo $TOWORKON | cut -f 5 -d"_" | cut -f 1 -d".")
      echo " I am working on date: $DATEON" 

      # create the outfile name
      TOOUT="${NEW_TPXOEXTR_PRE}_${DATEON}.nc"
      echo "New file ${TOOUT}"

      # mask the field
      python SeaOverLand.py ${TOWORKON} ${TPXOMASK} ${MESHMASK} ${TOOUT} ${LOOPNUM}
      #bsub -q s_long -n 1 -P 0284 python SeaOverLand.py ${TOWORKON} ${TPXOMASK} ${MESHMASK} ${TOOUT} ${LOOPNUM} 
      echo "Running $DATEON .."

   done 
done

echo "Done.."

###########################
