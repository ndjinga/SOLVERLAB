#!/bin/bash


if [ ${1:-none} == "-h" ]; then
  echo "
Manage and clean current PACKAGESPY directory to get
lighweighted PACKAGESPY for SOLVERLAB standalone installation
(i.e. without MATIX codes developments)

usage:
> cd .../PACKAGESPY    # mandatory
> ./configure_solverlabGUI.sh targetName

targetName: [ solverlab_standalone |
              solverlab_in_salome |
              solverlab_in_matix ]

example:
> ./configure_solverlabGUI.sh solverlab_in_salome
> ./configure_solverlabGUI.sh solverlab_in_matix
"
  exit 0  # ok
fi

# default parameter for script as environ variables 'conf_...' (as policy)
# may be more parameters later
export conf_target=${1:-solverlab_in_matix}  # keep maximum as default
# may be from sat later export conf_target=${conf_branch:-solverlab_in_matix}  # keep maximum as default

##### utilities #####

RED="\e[0;31m"
GREEN="\e[0;32m"
YELLOW="\e[0;33m"
BOLD="\e[1m"
NC="\e[0m" # No Color

function f_red {
  echo -e ${RED}${@}${NC}
}

function f_green {
  echo -e ${GREEN}${@}${NC}
}

function f_yellow {
  echo -e ${YELLOW}${@}${NC}
}

function f_error {
  f_red "\nERROR: "${@}"\n"
  exit 1  # ko
}

function f_warning {
  f_red "\nWARNING: "${@}"\n"
}

function f_delete_force {
  [ -z ${1} ] && f_error "f_delete_force without parameter, fit it!"
  [ -e ${1} ] || f_warning "inexisting file/directory '$1'"
  [ -d ${1} ] && (echo "rm -rf $1" ; rm -rf $1)
  [ -f ${1} ] && (echo "rm -f $1" ; rm -f $1)
}

function f_delete_dir {
  [ -d ${1:-none} ] || f_warning "inexisting directory $1" && f_delete_force $1
}

function f_delete_file {
  [ -f ${1:-none} ] || f_warning "inexisting file $1" && f_delete_force $1
}

function f_log_git {
  git remote -v | head -1 > git_history.tmp
  git --no-pager log -1 --graph --pretty=format:"%h %ad | %s%d <%an>" --graph --date=short >> git_history.tmp
  cat git_history.tmp
}

function f_log_result {	
  f_green '\nresult'
  tree -d -L 2
	f_warning "TODO some setenv.sh and else..."
}


###### main #####

# set -x   # echoes only for bash debug

#echo -e ${GREEN}
#env | grep -e '^conf_'  # beginning of line
#echo -e ${NC}
f_green "conf_target=${conf_target}"

# test if script launched in his own directory as current directory, is a bash posix headache
workdir=$(basename $(pwd))
[ ${workdir} == "PACKAGESPY_TULEAP" ] || f_error "current directory have to be '.../PACKAGESPY_TULEAP'"


if [ ${conf_target} == "solverlab_standalone" ]; then
	f_green '\nbegin clean dir(s) solverlab_standalone, keep the minimum'
	f_log_git
	f_delete_dir .git
	f_delete_dir OBSOLETES
	f_delete_dir pythonAppliMatix/amitexpy
	f_delete_dir pythonAppliMatix/cmdcpy
	f_delete_dir pythonAppliMatix/combspy
	f_delete_dir pythonAppliMatix/crescendopy
	f_delete_dir pythonAppliMatix/dartpy
	f_delete_dir pythonAppliMatix/doepy
	f_delete_dir pythonAppliMatix/ekinoxpy
	f_delete_dir pythonAppliMatix/microgenpy
	f_delete_dir pythonAppliMatix/neperpy
	f_delete_dir pythonAppliMatix/numodispy
	f_delete_dir pythonAppliMatix/octreepy
	f_delete_dir pythonAppliMatix/openglpy
	f_delete_dir pythonAppliMatix/ribbonpy
	f_delete_dir pythonAppliMatix/tesspy
	f_delete_dir pythonAppliMatix/uraniepy
	f_delete_dir pythonAppliMatix/voronoipy

	f_green "\nbegin clean file(s)"
	# f_delete_file .gitignore
	f_delete_file compile_PACKAGEPY_MATIX.py

	f_log_result


# almost idem solverlab_standalone for the moment
elif [ ${conf_target} == "solverlab_in_salome" ]; then
	f_green '\nbegin clean dir(s) solverlab_in_salome_light, keep the minimum'
	f_log_git
	f_delete_dir .git
	f_delete_dir OBSOLETES
	f_delete_dir pythonAppliMatix/amitexpy
	f_delete_dir pythonAppliMatix/cmdcpy
	f_delete_dir pythonAppliMatix/combspy
	f_delete_dir pythonAppliMatix/crescendopy
	f_delete_dir pythonAppliMatix/dartpy
	f_delete_dir pythonAppliMatix/doepy
	f_delete_dir pythonAppliMatix/ekinoxpy
	f_delete_dir pythonAppliMatix/microgenpy
	f_delete_dir pythonAppliMatix/neperpy
	f_delete_dir pythonAppliMatix/numodispy
	f_delete_dir pythonAppliMatix/octreepy
	f_delete_dir pythonAppliMatix/openglpy
	f_delete_dir pythonAppliMatix/ribbonpy
	f_delete_dir pythonAppliMatix/tesspy
	f_delete_dir pythonAppliMatix/uraniepy
	f_delete_dir pythonAppliMatix/voronoipy

	f_green "\nbegin clean file(s)"
	# f_delete_file .gitignore

	f_log_result


elif [ ${conf_target} == "solverlab_in_matix" ]; then
	f_green '\nbegin clean dir(s) solverlab_in_matix, keep the maximum'
	f_log_git
	echo f_delete_dir .git
	echo f_delete_dir OBSOLETES

	echo f_green '\nbegin clean file(s)'
	# echo f_delete_file .gitignore

	f_log_result

else
  f_error 'Unexpected $conf_target='${conf_target}

fi

echo
f_green ${0} seems to be OK

exit 0
