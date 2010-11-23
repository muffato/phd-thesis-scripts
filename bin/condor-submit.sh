#!/bin/bash

# Test nb arguments
if [ $# = 0 ]
then
	echo
	echo "Usage: $0 commande_classique avec arguments" 1>&2
	echo "Commandes speciales:  ! FICHIER_STDIN" 1>&2
	echo "                      / FICHIER_STDOUT" 1>&2
	echo "                     // FICHIER_STDERR" 1>&2
	echo "                    /// CONDOR_LOG" 1>&2
	echo "                     \~ REQUIREMENTS" 1>&2
	echo "                      % NB_ITER" 1>&2
	echo "                      + PRIORITY" 1>&2
	echo "                      $ [si NiceUser]" 1>&2
	echo "                      = [pour ne pas executer]" 1>&2
	echo "                     \* [pour utiliser le wrapper]" 1>&2
	echo '                      ^ [Notification = Always > Complete > Error > *Never*]' 1>&2
	echo "Exemple de REQUIREMENTS: '(machine != \"heimdall.ens.fr\") || (slotid <= 3)'"
	exit 1
fi

# Repertoire des fichiers de condor
CONDOR_WORK_DIR="$HOME/condor/"
BASENAME=$CONDOR_WORK_DIR`date +%F.%T.%N | sed 's/:/-/g'`

PROG=`which $1`
ARGS=
INPUT=
OUTPUT=$BASENAME".output"
LOG=$BASENAME".log"
CONDOR_LOG=$BASENAME".condor_log"
REQUIREMENTS='(machine != "heimdall.ens.fr")'
NICEUSER=False
NBQUEUE=
NOTIFICATION='never'
EXECUTE='YES'
PRIORITY=
WRAPPER=

shift
while [ $# -ge 1 ]
do
	if [ "$1" == "!" ]
	then
		INPUT=$2
		shift
	elif [ "$1" == "/" ]
	then
		OUTPUT=$2
		shift
	elif [ "$1" == "//" ]
	then
		LOG=$2
		shift
	elif [ "$1" == "///" ]
	then
		CONDOR_LOG=$2
		shift
	elif [ "$1" == "*" ]
	then
		WRAPPER="YES"
	elif [ "$1" == "~" ]
	then
		REQUIREMENTS=$2
		shift
	elif [ "$1" == "$" ]
	then
		NICEUSER=True
	elif [ "$1" == "%" ]
	then
		QUEUE=$2
		shift
	elif [ "$1" == "+" ]
	then
		PRIORITY=$2
		shift
	elif [ "$1" == "^" ]
	then
		NOTIFICATION=$2
		shift
	elif [ "$1" == "=" ]
	then
		EXECUTE='NO'
	else
		ARGS="$ARGS $1"
	fi
	shift
done

if [ "$WRAPPER" == "YES" ]
then
	ARGS="$PROG $ARGS"
	PROG="$HOME/bin/wrapper"
fi


##########################
#soumission du job condor#
##########################


if [ "$EXECUTE" == "YES" ]
then
	CONDOR_TASK_FILE=$BASENAME".condor_task"
else
	CONDOR_TASK_FILE=""
fi

echo "
# $CONDOR_TASK_FILE
Executable = $PROG
Universe = vanilla
Output = $OUTPUT
Input = $INPUT
Error = $LOG
Log = $CONDOR_LOG
Arguments = $ARGS
GetEnv = True
Initialdir = `pwd`
should_transfer_files = NO
Requirements = $REQUIREMENTS
Notification = $NOTIFICATION
NiceUser = $NICEUSER
Priority = $PRIORITY
Rank = kflops+1000*Memory
Queue $QUEUE
" | tee $CONDOR_TASK_FILE

if [ "$EXECUTE" == "YES" ]
then
	condor_submit $CONDOR_TASK_FILE 1>&2
fi

