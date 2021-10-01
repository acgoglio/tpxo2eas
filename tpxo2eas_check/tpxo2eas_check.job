#!/bin/bash
#
#BSUB -J tpxo_check          # Name of the job.
#BSUB -o /work/oda/ag15419/job_scratch/extr_%J.out  # Appends std output to file %J.out.
#BSUB -e /2ork/oda/ag15419/job_scratch/extr_%J.err  # Appends std error to file %J.err.
#BSUB -q s_medium
#BSUB -n 1    # Number of CPUs
#BSUB -P 0284
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 02/03/2021
#
#set -u
set -e
#set -x 
########################
echo "*********** TPXO extraction check  *********"

module load anaconda/3.7 curl/7.70.0 cmake/3.17.3 gams/28.2.0 gcc_9.1.0/9.1.0 gcc_9.1.0/gempack/12.885 gcc_9.1.0/OpenBLAS/0.3.9 gcc_9.1.0/papi/6.0.0 gcc_9.1.0/R/3.6.1 modules mysql/5.7.28 ncl/6.6.2 sqlite/3.32.2 subversion/1.14.0 wgrib/1.8.1.0b impi20.1/19.7.217 impi20.1/esmf/8.0.1-intelmpi-64-g impi20.1/hdf5/1.12.0 impi20.1/hdf5-threadsafe/1.12.0 impi20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/papi/6.0.0 impi20.1/parallel-netcdf/1.12.1 impi20.1/petsc/3.13.2 impi20.1/zoltan/3.8 intel20.1/20.1.217 intel20.1/advisor intel20.1/boost/1.73.0 intel20.1/cdo/1.9.8 intel20.1/cnvgrib/3.1.1 intel20.1/eccodes/2.17.0 intel20.1/esmf/8.0.1-mpiuni-64-g intel20.1/esmf/8.0.1-mpiuni-64-O intel20.1/exactextract/545f0d6 intel20.1/g2lib/3.1.0 intel20.1/gdal/3.1.0 intel20.1/hdf5/1.12.0 intel20.1/hdf5-threadsafe/1.12.0 intel20.1/inspector intel20.1/itac intel20.1/libemos/4.5.9 intel20.1/libemos/4.5.9 intel20.1/magics/3.3.1 intel20.1/nco/4.9.3 intel20.1/ncview/2.1.8 intel20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/proj/7.0.1 intel20.1/R/4.0.2 intel20.1/szip/2.1.1 intel20.1/udunits/2.2.26 intel20.1/valgrind/3.16.0 intel20.1/vtune intel20.1/w3lib/2.0.6 intel20.1/wgrib2/2.0.8


WORKDIR="/work/oda/ag15419/tmp/new_tpxoextr/"
TPXODB_PATH="/data/oda/ag15419/TPXO9_DATA/new_tpxoextr/" #"/data/opa/mfs-dev/Med_static/MFS_TPXO_V1/"
POINT_LAT=41.895832
POINT_LON=16.166666

cd ${WORKDIR}

for EXTRYEAR in 2022 ; do # 2019 2020 2021 2022 2023; do
  IDX_NC=0
  for MONTHS_IDX in 01 02 03 04 05 06 07 08 09 10 11 12; do
    # Check min, mean and max on the whole domain
    echo "# TPXO EXTRACTION YEAR ${EXTRYEAR} MONTH ${MONTHS_IDX}"
    echo "# TPXO EXTRACTION YEAR ${EXTRYEAR} MONTH ${MONTHS_IDX}" > ${WORKDIR}/infon_${EXTRYEAR}${MONTHS_IDX}.txt
    for TSEXT in $(ls ${TPXODB_PATH}/${EXTRYEAR}/tpxoextr_${EXTRYEAR}${MONTHS_IDX}??.nc) ; do
        cdo infon $TSEXT | grep "tide_z" >> ${WORKDIR}/infon_${EXTRYEAR}${MONTHS_IDX}.txt
    done
    GPL_FILE=${EXTRYEAR}${MONTHS_IDX}_plot.gpl
    cat << EOF > ${WORKDIR}/${GPL_FILE}
set term pngcairo size 1700,700 font "Times-New-Roman,16"
set output "tpxo_mmm_${EXTRYEAR}${MONTHS_IDX}.png"
set multiplot layout 3,1 title "TPXO extracted data YEAR ${EXTRYEAR} MONTH ${MONTHS_IDX}"
#set title "TPXO extracted data YEAR ${EXTRYEAR} MONTH ${MONTHS_IDX}"
set xlabel "Date"
set xdata time
set timefmt "%Y-%m-%d %H:%M:%S"
#set xrange ["${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-01":"${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}"]
set format x "%d/%m/%Y"
set ylabel "Sea Level [m]"
set grid
set key Left
set datafile missing "nan"
plot 'infon_${EXTRYEAR}${MONTHS_IDX}.txt' using 3:9 with line lw 3 lt rgb '#1f77b4' title "Domain TPXO min"
plot 'infon_${EXTRYEAR}${MONTHS_IDX}.txt' using 3:10 with line lw 3 lt rgb '#ff7f0e' title "Domain TPXO mean"
plot 'infon_${EXTRYEAR}${MONTHS_IDX}.txt' using 3:11 with line lw 3 lt rgb '#d62728' title "Domain TPXO max"
EOF
    gnuplot < ${WORKDIR}/${GPL_FILE} || echo "Problem with plot.. Why?"
 
    # Check time serie in single grid point
    GPL_FILE=${EXTRYEAR}${MONTHS_IDX}_pointplot.gpl
    cat << EOF > ${WORKDIR}/${GPL_FILE}
set term pngcairo size 1700,700 font "Times-New-Roman,16"
set output "tpxo_point_${EXTRYEAR}${MONTHS_IDX}.png"
set title "TPXO extracted data YEAR ${EXTRYEAR} MONTH ${MONTHS_IDX} in Lat/Lon ${POINT_LAT}/${POINT_LON}"
set xlabel "time-steps [2 min]"
#set xdata time
#set timefmt "%Y-%m-%d %H:%M:%S"
#set xrange ["${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-01":"${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}"]
#set format x "%d/%m/%Y"
set ylabel "Sea Level [m]"
set grid
set key Left
set datafile missing "nan"
plot 'tpxo2check_${EXTRYEAR}${MONTHS_IDX}.txt' using 0:1 with line lw 3 lt rgb '#2ca02c' title "TPXO in Lat/Lon ${POINT_LAT}/${POINT_LON}"
EOF
    gnuplot < ${WORKDIR}/${GPL_FILE} || echo "Problem with plot.. Why?"
  done
done

