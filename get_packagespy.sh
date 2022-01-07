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
get_branch=${get_branch:-ym_M30}
get_configure=${get_configure:-configure_solverlabGUI.sh}

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

envs | grep "get_"

# test if script launched in his own directory as current directory
workdir=$(basename $(pwd))
[ ${workdir} == "SOLVERLAB" ] || f_error "current directory have to be '.../SOLVERLAB'"

# clone packagespy under .../SOLVERLAB, packagespy have not to already exist
[ -d "packagespy" ] && f_error "directory .../SOLVERLAB/packagespy existing yet, fix it (remove ?)."
git clone --branch ${get_branch} https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git || f_error "git clone packagespy"

# this local copy is only for EZ development of packagespy/${get_configure}
[ ${USER:-none} == "wambeke" ] && cp ${get_configure} packagespy/.

# store some current packagespy git info before configure/clean
cd packagespy && f_git_log | tee git_history.txt

# theorically something as 'packagespy/configure_solverlabGUI.sh'
[ -f "$get_configure" ] || f_error "inexisting file packagespy/$get_configure"
${get_configure} || f_error "problem in packagespy/$get_configure"

echo
f_green ${0} seems to be OK

exit 0
