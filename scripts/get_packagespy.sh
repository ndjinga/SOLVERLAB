#!/bin/bash

if [ ${1:-none} == "-h" ]; then
  echo "usage:

> cd .../SOLVERLAB
> [parameters] get_packagespy.sh

set default parameter(s) as environ variable(s), if needed
example:
> get_branch=master get_packagespy.sh
"
  exit 0
fi

##### default parameter(s) for script as environ variables 'get_...' (as policy)
#  may be more later
export get_base=${get_base:-"https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git"}
export get_branch=${get_branch:-ym_M30}
export get_configure=${get_configure:-configure_solverlabGUI.sh}
export get_mode_devel=${get_mode_devel:-OFF}  # if ON then keep PACKAGESPY with git (for developer)

##### utilities #####

RED="\e[0;31m"
GREEN="\e[0;32m"
BOLD="\e[1m"
NC="\e[0m" # No Color

function f_red {
  echo -e ${RED}${@}${NC}
}

function f_green {
  echo -e ${GREEN}${@}${NC}
}

function f_error {
  f_red "\nERROR: "${@}"\n"
  exit 1
}

function f_git_log {
  # see small commits/tag history
  echo -e "\ngit log:"
  git --no-pager log -5 --graph --pretty=format:"%h %ad | %s%d <%an>" --date=short
  echo -e "\n\ngit branches:"
  git branch -a
  echo -e "\ngit status:"
  git status
}


###### main #####

# echoes only for debug
# set -x

f_green "parameters:"
env | grep -e '^get_' ; echo  # beginning of line

# test if script launched in his own directory as current directory
workdir=$(basename $(pwd))
[ ${workdir} == "SOLVERLAB" ] || f_error "current directory have to be '.../SOLVERLAB'"

# clone packagespy.git under .../SOLVERLAB, PACKAGESPY have not to already exist
[ -d "PACKAGESPY" ] && f_error "directory .../SOLVERLAB/PACKAGESPY existing yet, fix it (remove ?)."
git clone --branch ${get_branch} ${get_base} PACKAGESPY || f_error "git clone PACKAGESPY"

# this local copy is only for EZ development of packagespy/${get_configure}
# [ ${USER:-none} == "wambeke" ] && cp ${get_configure} PACKAGESPY/.

# store some current packagespy git info before configure/clean
cd PACKAGESPY && f_git_log | tee git_history.txt

# theorically something as 'PACKAGESPY/configure_solverlabGUI.sh'
[ -f "$get_configure" ] || f_error "inexisting file PACKAGESPY/${get_configure}"

if [ ${get_mode_devel} == OFF ]; then
  ${get_configure} || f_error "problem in PACKAGESPY/$get_configure"
else
  f_red "\nget_mode_devel=${get_mode_devel}  # keep PACKAGESPY .git for developments"
fi

echo
f_green ${0} seems to be OK

exit 0
