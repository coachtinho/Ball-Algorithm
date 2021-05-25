#!/bin/bash

# This script can be run anywhere if you set the variable below

# It has the following options:
# no-color: print output uncolored
# no-compact: print the output of each test to the console
# no-clean: keep all the logs

# Set the path to the folder with the tests
TESTS_PATH="tests"
# Set the path to the folder where the logs are going to be stored
LOG_FOLDER="log"
# Set the path to the program from the tests folder
PROG="../ballAlg-mpi"
QUERY="../ballQuery"
#
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# DO NOT TOUCH AFTER THIS LINE
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

CLEAN=true
COMPACT=false
COLOR=true

if [ $# -ne 0 ]; then
    for i in $@
    do
        case ${i,,} in
            "no-clean")
                CLEAN=false
                ;;
            "no-compact")
                COMPACT=false
                ;;
            "no-color")
                COLOR=false
                ;;
            *)
                echo "Unknown option ${BOLD}$i${RESET}"
                ;;
        esac
    done
fi


cd $TESTS_PATH

if [ "$COLOR" = true ]; then
    RESET="\e[0m"
    BOLD="\e[1m"
    BLUE="\e[34m"
    ORANGE="\e[33m"
    RED="\e[31m"
    GREEN="\e[32m"
fi

if [ ! -d $LOG_FOLDER ]; then
    mkdir $LOG_FOLDER
fi

function clean_lines_up {
    for i in $(seq $1)
    do
        echo -e -n "\033[A\r\033[K"
    done
}

SUCCESS=0
TOTAL=$(ls *.in | wc -l)
CURRENT=0
for file in *.in
do
    CURRENT=$(($CURRENT + 1))
    util=${file%.in}
    echo "FILE: ${file} (${CURRENT}/${TOTAL})"
    echo -n "Executing program... "
    rm expected/${util}.query.mine &> /dev/null

    if [ $(echo $PROG | grep mpi) ]; then
        ${QUERY} <(mpirun --use-hwthread-cpus -n 4 $PROG $(cat ${file}) 2> /dev/null) $(cat ${util}.query) 2>/dev/null > expected/${util}.query.mine
    else
        ${QUERY} <($PROG $(cat ${file}) 2> /dev/null) $(cat ${util}.query) 2>/dev/null > expected/${util}.query.mine
    fi

    if [ $? -eq 0 ]; then
        echo "DONE"
    else
        echo -e "FAILED\n"
        if [ "$COMPACT" = true ]; then
            clean_lines_up 2
        fi
        continue
    fi

    echo -n "Comparing outputs... "
    if [ "$(diff -q -b expected/${util}.query.out expected/${util}.query.mine)" != "" ]; then
        echo -e "DIFFERENT\n"
        if [ "$COMPACT" = true ]; then
            clean_lines_up 3
        fi
        continue
    else
        echo "EQUAL"
    fi


    echo -e "Cleaning...\n "
    rm expected/${util}.query.mine

    SUCCESS=$(($SUCCESS + 1))
    if [ "$COMPACT" = true ]; then
        clean_lines_up 4
    fi
done
echo -e "${BOLD}${BLUE}FINAL RESULTS:${RESET}"
echo -e "${BOLD}Passed: ${GREEN}${SUCCESS}${RESET}"
echo -e "${BOLD}Failed: ${RED}$(($TOTAL - $SUCCESS))${RESET}"

perc=$(echo "scale=2; $SUCCESS * 100 / $TOTAL" | bc -l)

if [ $(echo "$perc > 25" | bc -l) -eq 0 ]; then
    color=$RED
elif [ $(echo "$perc > 50" | bc -l) -eq 0 ]; then
    color=$YELLOW
elif [ $(echo "$perc > 80" | bc -l) -eq 0 ]; then
    color=$ORANGE
else
    color=$GREEN
fi

echo -e "${BOLD}Percentage: ${color}${perc}% ($SUCCESS / $TOTAL)${RESET}"

if [ $SUCCESS -ne $TOTAL ]; then
    # echo -e "${BOLD}${ORANGE}NOTE:${RESET} Check ${BOLD}${TESTS_PATH}/${LOG_FOLDER}${RESET} to see the full results"
    exit 1
fi
