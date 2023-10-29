#!/bin/bash

# This script gets the strings from code to internationalize from the source code

I18HOME=`dirname $0`
SRC_DIR=$I18HOME/../..

# get strings for french translation
echo "Build strings for French, create and merging solverlabGUI.po"

poFile=$I18HOME/fr/LC_MESSAGES/solverlabGUI.po
refFile=$I18HOME/fr/LC_MESSAGES/ref.pot

cp ${poFile} ${poFile}_old

xgettext $SRC_DIR/src/i18n/*.py \
         $SRC_DIR/*.py \
         $SRC_DIR/commands/*.py \
         $SRC_DIR/solverlabpy/*.py \
         --no-wrap \
         --no-location \
         --language=Python \
         --omit-header \
         --output=${refFile}

msgmerge --quiet --update --previous $poFile $refFile

#retirer les messages obsoletes « #~ »
#msgattrib --no-obsolete -o $poFile $poFile

#ne pas retirer les messages obsolètes « #~ »
msgattrib --previous --output-file ${poFile} ${poFile}

rm $refFile

echo "Do translate stuff please..."
meld ${poFile} ${poFile}_old

echo "Do not forget 'translate.py' or 'translate.sh' to create solverlabGUI.mo"

