#!/bin/bash
set -e
#set -x 
################### ENV SETTINGS ##############################
SRC_DIR="/users_home/oda/ag15419/tpxo2eas/"
echo "Job-task element: ${LSF_PM_TASKID}"
################### PREPROC ###################################

# Source ini file
  INIFILE="${SRC_DIR}/tpxo_dayextr.ini"
  source $INIFILE
  echo "source $INIFILE ... Done!"

  module load ${EXTR_MODULE}  
  #module list

# Read and check infos (work dir, file names, archive dir, etc.)

  # Workdir check and subdir mk
  if [[ -d $EXTR_WORKDIR ]]; then
     cd $EXTR_WORKDIR
     #mkdir day_${LSF_PM_TASKID}
     #cd day_${LSF_PM_TASKID}
     #EXTR_SUBWORKDIR=${EXTR_WORKDIR}/day_${LSF_PM_TASKID}
     echo "WORKDIR: $(pwd)"
     
     # Clean workdir
     #echo "WARNING: I am going to remove all files in $EXTR_WORKDIR ..."
     #sleep 10
     #for TO_BE_RM in $( ls $EXTR_WORKDIR ); do
     ##    rm -vr $EXTR_WORKDIR/$TO_BE_RM
     #    echo $TO_BE_RM
     #done
     # Cp the exe file to the workdir
     if [[ ! -f ${EXE_OTPS} ]]; then 
        echo "I need to copy the exe file to the workdir.."
        cp ${SRC_DIR}/${EXE_OTPS} .
     fi

  else
     echo "ERROR: WORKDIR $EXTR_WORKDIR NOT FOUND!!"
     exit
  fi

  # Cp the Model_atlas file to the workdir
  if [[ ! -d DATA ]] && [[ ! -f DATA/Model_atlas ]]; then
        echo "I need to copy the DATA dir to the workdir.."
        cp -rv ${SRC_DIR}/DATA .
  #else
   #     echo "ERROR: DATA/Model_atlas NOT FOUND in SRC_DIR.."
  fi



  # Set the date based on the index of the task
  SET_OF_TASK=$(( ${1} -1 ))
  SET_OF_TASK=$(( ${SET_OF_TASK} *36 ))
  echo "set_of_task: $SET_OF_TASK"
  DAY_IDX=$(( ${LSF_PM_TASKID} + ${SET_OF_TASK} ))
  echo "The extraction starts on date: ${EXTR_STARTDATE}"
  DAY2EXTR=$(date -d "${EXTR_STARTDATE} ${DAY_IDX} day" +%Y%m%d)  
  echo "This task element extracts date: ${DAY2EXTR}"

  # Check mesh_mask file
  MESHMASK="${MESHMASK_PATH}/${MESHMASK_FILE}"
  if [[ -f ${MESHMASK} ]]; then
     echo "Lat and lon are taken from file: ${MESHMASK}"
  else
     echo "ERROR: mesh mask file NOT found..Why?"
     exit
  fi

  # Check out dir and built outfile name
  if [[ -d ${OUTNC_PATH} ]]; then
     OUTNC_NAME=$( echo "$OUTNC_TEMPLATE" | sed -e "s/%YYYYMMDD%/${DAY2EXTR}/g" )
     OUTFILE_NC="${OUTNC_PATH}/${OUTNC_NAME}"
     echo "The output.nc file is: ${OUTFILE_NC}"
  else
     echo "ERROR: OUT dir NOT found..Why?"
     exit
  fi

  # Run daily extraction..
  if [[ -f ${EXE_OTPS} ]]; then 
     echo "Run the script with the following args: "
     echo "./${EXE_OTPS} ${DAY2EXTR} ${MESHMASK} ${OUTFILE_NC}"
     time ./${EXE_OTPS} ${DAY2EXTR} ${MESHMASK} ${OUTFILE_NC}
  else
     echo "ERROR: Not found exe file ${EXE_OTPS}...Why?!"
     exit
  fi


###################### POSTPROC ###########################

# Output check

# Archive


