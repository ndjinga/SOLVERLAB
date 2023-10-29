#!/bin/bash

# EZ install splash matix in INSTALL/SALOME/share/salome/resources/salome
# to avoid manage useless MATIX_PROFILE (etc) module from MATIX 3.0 211201
# and use directly SALOME module

if [ -z "$PACKAGESPY_ROOT_DIR" ]
then
  echo "You need set environ PACKAGESPY_ROOT_DIR before (...as matix context)."
  exit 1
fi

set -x

cp ${PACKAGESPY_ROOT_DIR}/splash/icon_splash_matix.png  ${PACKAGESPY_ROOT_DIR}/../SALOME/share/salome/resources/salome/splash.png
cp ${PACKAGESPY_ROOT_DIR}/splash/icon_applogo_matix.png ${PACKAGESPY_ROOT_DIR}/../SALOME/share/salome/resources/salome/app_logo.png

