#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
highlight -h
USAGE: highlight [OPTIONS]... [FILES]...

General options:

 -B, --batch-recursive=<wc>     convert all matching files, searches subdirs
                                  (Example: -B '*.cpp')
 -D, --data-dir=<directory>     set path to data directory
     --add-config-dir=<dir>     set path to an additional config directory
     --config-file=<file>       set path to a lang or theme file
 -d, --outdir=<directory>       name of output directory
 -h, --help                     print this help
 -i, --input=<file>             name of single input file
 -o, --output=<file>            name of single output file
 -p, --list-langs               list installed language definitions
 -P, --progress                 print progress bar in batch mode
 -q, --quiet                    supress progress info in batch mode
 -S, --syntax=<type>            specify type of source code
 -v, --verbose                  print debug info
 -w, --list-themes              list installed colour themes
     --force                    generate output if input syntax is unknown
     --plug-in=<script>         execute Lua plug-in script; repeat option to
                                  execute multiple plug-ins
     --plug-in-read=<path>      set input file for a plug-in (e.g. "tags")
     --print-config             print path configuration
     --print-style              print stylesheet only (see --style-outfile)
     --skip=<list>              ignore listed unknown file types
                                  (Example: --skip='bak;c~;h~')
     --start-nested=<lang>      define nested language which starts input
                                  without opening delimiter
     --validate-input           test if input is text, remove Unicode BOM
     --version                  print version and copyright information


Output formatting options:

 -O, --out-format=<format>      output file in given format
                                  <format>=[html, xhtml, latex, tex,
                                  odt, rtf, ansi, xterm256, bbcode, svg]
 -c, --style-outfile=<file>     name of style file or print to stdout, if
                                  'stdout' is given as file argument
 -e, --style-infile=<file>      file to be included in style-outfile
 -f, --fragment                 omit document header and footer
 -F, --reformat=<style>         reformats and indents output in given style
                                  <style>=[allman, banner, gnu,
                                  horstmann, java, kr, linux, otbs,
                                  stroustrup, whitesmith]
 -I, --include-style            include style definition
 -J, --line-length=<num>        line length before wrapping (see -W, -V)
 -j, --line-number-length=<num> line number width incl. left padding
 -k, --font=<font>              set font (specific to output format)
 -K, --font-size=<num?>         set font size (specific to output format)
 -l, --line-numbers             print line numbers in output file
 -m, --line-number-start=<cnt>  start line numbering with cnt (assumes -l)
 -s, --style=<style>            set colour style (see -w)
 -t, --replace-tabs=<num>       replace tabs by <num> spaces
 -T, --doc-title=<title>        document title
 -u, --encoding=<enc>           set output encoding which matches input file
                                  encoding; omit encoding info if set to NONE
 -V, --wrap-simple              wrap long lines without indenting function
                                  parameters and statements
 -W, --wrap                     wrap long lines
     --wrap-no-numbers          omit line numbers of wrapped lines
                                  (assumes -l)
 -z, --zeroes                   pad line numbers with 0's
     --kw-case=<case>           change case of case insensitive keywords
                                  <case> =  [upper, lower, capitalize]
     --delim-cr                 set CR as end-of-line delimiter (MacOS 9)
     --no-trailing-nl           omit trailing newline


(X)HTML output options:

 -a, --anchors                  attach anchor to line numbers
 -y, --anchor-prefix=<str>      set anchor name prefix
 -N, --anchor-filename          use input file name as anchor prefix
 -C, --print-index              print index with hyperlinks to output files
 -n, --ordered-list             print lines as ordered list items
     --class-name=<name>        set CSS class name prefix;
                                  omit class name if set to NONE
     --inline-css               output CSS within each tag (verbose output)
     --enclose-pre              enclose fragmented output with pre tag 
                                  (assumes -f)


LaTeX output options:

 -b, --babel                    disable Babel package shorthands
 -r, --replace-quotes           replace double quotes by \dq{}
     --pretty-symbols           improve appearance of brackets and other symbols


RTF output options:

 -x, --page-size=<ps>           set page size 
                                  <ps> = [a3, a4, a5, b4, b5, b6, letter]
     --char-styles              include character stylesheets


SVG output options:

     --height                   set image height (units allowed)
     --width                    set image width (see --height)


GNU source-highlight compatibility options:

     --doc                      create stand alone document
     --no-doc                   cancel the --doc option
     --css=filename             the external style sheet filename
     --src-lang=STRING          source language
 -t, --tab=INT                  specify tab length
 -n, --line-number[=0]          number all output lines, optional padding
     --line-number-ref[=p]      number all output lines and generate an anchor,
                                  made of the specified prefix p + the line
                                  number  (default='line')
     --output-dir=path          output directory
     --failsafe                 if no language definition is found for the
                                  input, it is simply copied to the output


If no in- or output files are specified, stdin and stdout will be used.
HTML will be generated unless an other output format is given. Style definitions
are stored in highlight.css (HTML, XHTML, SVG) or highlight.sty (LaTeX, TeX)
if neither -c nor -I is given.
Reformatting code (-F) will only work with C, C++, C# and Java input files.
Wrapping lines with -V or -W will cause faulty highlighting of long single
line comments and directives. Use with caution.

Updates and information: http://www.andre-simon.de/
"""

"""
highlight --list-themes

Installed themes (located in /usr/share/highlight/themes/):

acid
aiseered
andes
anotherdark
autumn
baycomb
bclear
biogoo
bipolar
blacknblue
bluegreen
breeze
bright
camo
candy
clarity
dante
darkblue
darkbone
darkness
darkslategray
darkspectrum
denim
dusk
earendel
easter
edit-anjuta
edit-eclipse
edit-emacs
edit-flashdevelop
edit-gedit
edit-jedit
edit-kwrite
edit-matlab
edit-msvs2008
edit-nedit
edit-vim-dark
edit-vim
edit-xcode
ekvoli
fine_blue
freya
fruit
golden
greenlcd
kellys
leo
lucretia
manxome
maroloccio
matrix
moe
molokai
moria
navajo-night
navy
neon
night
nightshimmer
nuvola
olive
orion
oxygenated
pablo
peaksea
print
rand01
rdark
relaxedgreen
rootwater
seashell
solarized-dark
solarized-light
tabula
tcsoft
vampire
whitengrey
xoria256
zellner
zenburn
zmrok

Use name of the desired theme with the --style option.
"""

"""
highlight --list-langs

Installed language definitions (located in /usr/share/highlight/langDefs/):

ABAP/4                        : abap4 ( abp )
ABC                           : abc
Advanced Backus-Naur Form     : abnf
Action Script                 : actionscript ( as )
ADA95                         : ada ( a adb ads gnad )
Agda                          : agda
ALGOL 68                      : algol ( alg )
AMPL                          : ampl ( dat run )
AMTrix                        : amtrix ( hnd s4 s4h s4t t4 )
AppleScript                   : applescript
Arc                           : arc
ARM                           : arm
AS/400 CL                     : as400cl
ASCEND                        : ascend ( a4c )
ASP                           : asp ( asa )
Abstract                      : aspect ( was wud )
Assembler                     : assembler ( asm )
Applied Type System           : ats ( dats )
AutoHotKey                    : autohotkey ( ahk )
AutoIt                        : autoit ( au3 )
Avenue                        : avenue
(G)AWK                        : awk
DOS Batch                     : bat ( cmd )
BBcode                        : bbcode
BCPL                          : bcpl
BibTeX                        : bibtex ( bib )
Biferno                       : biferno ( bfr )
Bison                         : bison ( y )
Blitz Basic                   : blitzbasic ( bb )
BM Script                     : bms
Backus-Naur Form              : bnf
Boo                           : boo
C and C++                     : c ( c++ cc cpp cu cxx h hh hpp hxx )
Ceylon                        : ceylon
Charmm                        : charmm ( inp )
CHILL                         : chill ( chl )
Clean                         : clean ( icl )
ClearBasic                    : clearbasic ( cb )
Clipper                       : clipper
Clojure                       : clojure
Clips                         : clp
COBOL                         : cobol ( cbl cob )
ColdFusion MX                 : coldfusion ( cfc cfm )
Crack                         : crk
C#                            : csharp ( cs )
CSS                           : css
D                             : d
Dart                          : dart
Diff                          : diff ( patch )
Dylan                         : dylan
Extended Backus-Naur Form     : ebnf
Eiffel                        : eiffel ( e se )
Erlang                        : erlang ( erl hrl )
Euphoria                      : euphoria ( eu ew ex exw wxu )
Express                       : express ( exp )
FAME                          : fame ( fame )
Felix                         : felix ( flx )
Fortran 77                    : fortran77 ( f for ftn )
Fortran 90                    : fortran90 ( f90 f95 )
Frink                         : frink
F#                            : fsharp ( fs fsx )
Java FX                       : fx
Gambas                        : gambas ( class )
Go                            : go
Graphviz                      : graphviz ( dot )
Haskell                       : haskell ( hs )
haXe                          : haxe ( hx )
Hecl                          : hcl
HTML                          : html ( htm xhtml )
Apache Config                 : httpd
Icon                          : icon ( icn )
IDL                           : idl
Interactive Data Language     : idlang
Lua (for LuaTeX)              : inc_luatex
Informix                      : informix ( 4gl )
INI                           : ini
Inno Setup                    : innosetup ( iss )
INTERLIS                      : interlis ( ili )
IO                            : io
Jasmin                        : jasmin ( j )
Java                          : java ( groovy grv )
Javascript                    : js
JSP                           : jsp
LDAP                          : ldif
Haskell LHS                   : lhs
Lilypond                      : lilypond ( ly )
Limbo                         : limbo ( b )
Linden Script                 : lindenscript ( lsl )
Lisp                          : lisp ( cl clisp el lsp sbcl scom )
Logtalk                       : logtalk ( lgt )
Lotos                         : lotos
Lotus                         : lotus ( ls )
Lua                           : lua
Luban                         : luban ( lbn )
Make                          : make ( mak mk )
Maple                         : maple ( mpl )
Matlab                        : matlab ( m )
Maya                          : maya ( mel )
Mercury                       : mercury
Miranda                       : miranda
Modula2                       : mod2 ( def mod )
Modula3                       : mod3 ( i3 m3 )
Modelica                      : modelica ( mo )
MoonScript                    : moon
MaxScript                     : ms
MSSQL                         : mssql
Magic eXtensible Markup Language: mxml
Notation3 (N3), N-Triples, Turtle, SPARQL: n3 ( nt ttl )
Nasal                         : nasal ( nas )
NeXT Byte Codes               : nbc
Nemerle                       : nemerle ( n )
NetRexx                       : netrexx ( nrx )
Nice                          : nice
NSIS                          : nsis ( nsi )
Not eXactly C                 : nxc
Oberon                        : oberon ( ooc )
Objective C                   : objc
Objective Caml                : ocaml ( ml mli )
Octave                        : octave
Open Object Rexx              : oorexx
Object Script                 : os
Oz                            : oz
Paradox                       : paradox ( sc )
Pascal                        : pas
Perl                          : perl ( cgi perl pl plex plx pm )
PHP                           : php ( php3 php4 php5 php6 )
Pike                          : pike ( pmod )
PL/1                          : pl1 ( bdy ff fp fpp rpp sf sp spb spe spp sps wf wp wpb wpp wps )
PL/Perl                       : plperl
PL/Python                     : plpython
PL/Tcl                        : pltcl
POV-Ray                       : pov
Prolog                        : pro
Progress                      : progress ( i p w )
PostScript                    : ps
Microsoft PowerShell          : ps1
PATROL                        : psl
Pure                          : pure
Pyrex                         : pyrex ( pyx )
Python                        : python ( py )
Qore                          : q
QMake Project                 : qmake
Qu                            : qu
R                             : r
Rebol                         : rebol
Rexx                          : rexx ( rex rx the )
Relax NG                      : rnc
RPG                           : rpg
RPL Programming Language      : rpl
Ruby                          : ruby ( pp rb rjs ruby )
PowerPC Assembler             : s
SAS                           : sas
Scala                         : scala
Scilab                        : scilab ( sce sci )
Bash                          : sh ( bash ebuild eclass )
SMALL                         : small ( sma )
Smalltalk                     : smalltalk ( gst sq st )
Standard ML                   : sml
SNMP                          : snmp ( mib smi )
SNOBOL                        : snobol ( sno )
RPM Spec                      : spec
SPIN SQL                      : spn
PL/SQL                        : sql
Squirrel                      : squirrel ( nut )
Sybase SQL                    : sybase
Tcl/Tk                        : tcl ( itcl wish )
TCSH                          : tcsh
TeX and LaTeX                 : tex ( cls sty )
TypeScript                    : ts
Transact-SQL                  : tsql
TTCN3                         : ttcn3
Plain text                    : txt ( text )
UPC (and C, technically)      : upc
Vala                          : vala
Visual Basic                  : vb ( bas basic bi vbs )
Verilog                       : verilog ( v )
VHDL                          : vhd
XML                           : xml ( dtd ecf ent hdr hub jnlp nrm resx sgm sgml svg tld vxml wml xsd xsl )
SuperX++                      : xpp
Yaiff                         : yaiff
Zonnon                        : znn

Use name of the desired language with the --syntax option.
"""

"""
cat /etc/highlight/filetypes.conf

-- File extension and shebang mapping

FileMapping = {

 { Lang="ada",  Extensions={"adb", "ads", "a", "gnad"} },
 { Lang="algol",  Extensions={"alg"} },
 { Lang="ampl", Extensions={"dat", "run"} },
 { Lang="amtrix", Extensions={"s4", "s4t", "s4h", "hnd", "t4"} },
 { Lang="asm", Extensions={"a51", "29k", "68s", "68x", "x86"} },
 { Lang="asp", Extensions={"asa"} },
 { Lang="ats", Extensions={"dats"} },
 { Lang="aspect", Extensions={"was", "wud"} },
 { Lang="bat", Extensions={"cmd"} },
 { Lang="c", Extensions={"c++", "cpp", "cxx", "cc", "h", "hh", "hxx", "hpp", "cu"} },
 { Lang="charmm", Extensions={"inp"} },
 { Lang="coldfusion", Extensions={"cfc","cfm"} },
 { Lang="cobol", Extensions={"cob", "cbl"} },
 { Lang="diff", Extensions={"patch"} },
 { Lang="eiffel", Extensions={"e", "se"} },
 { Lang="erlang", Extensions={"hrl", "erl"} },
 { Lang="euphoria", Extensions={"ex", "exw", "wxu", "ew", "eu"} },
 { Lang="fortran77", Extensions={"f", "for", "ftn"} },
 { Lang="fortran90", Extensions={"f95", "f90"} },
 { Lang="gambas", Extensions={"class"} },
 { Lang="haskell", Extensions={"hs"} },
 { Lang="java", Extensions={"groovy", "grv"} },
 { Lang="limbo", Extensions={"b"} },
 { Lang="lisp", Extensions={"cl", "clisp", "el", "lsp", "sbcl", "scom"} },
 { Lang="make", Extensions={"mak", "mk"} },
 { Lang="snmp", Extensions={"mib", "smi"} },
 { Lang="ocaml", Extensions={"ml","mli"} },
 { Lang="mod2", Extensions={"mod", "def"} },
 { Lang="mod3", Extensions={"m3", "i3"} },
 { Lang="oberon", Extensions={"ooc"} },
 { Lang="php", Extensions={"php3", "php4", "php5", "php6"} },
 { Lang="pike", Extensions={"pmod"} },
 { Lang="pl1", Extensions={"ff", "fp", "fpp", "rpp","sf", "sp", "spb",
               "spp","sps", "wp", "wf", "wpp","wps","wpb","bdy","spe"} },
 { Lang="perl", Extensions={"pl","perl", "cgi", "pm", "plx", "plex"} },
 { Lang="progress", Extensions={"p", "i", "w"} },
 { Lang="ruby", Extensions={"rb","ruby", "pp", "rjs"} },
 { Lang="rexx", Extensions={"rex", "rx", "the"} },
 { Lang="sh", Extensions={"bash", "ebuild", "eclass"} },
 { Lang="smalltalk", Extensions={"st", "gst", "sq"} },
 { Lang="sybase", Extensions={"sp"} },
 { Lang="tcl", Extensions={"wish", "itcl"} },
 { Lang="tex", Extensions={"sty", "cls"} },
 { Lang="vb", Extensions={"bas", "basic", "bi", "vbs"} },
 { Lang="verilog", Extensions={"v"} },
 { Lang="html", Extensions={"htm", "xhtml"} },
 { Lang="xml", Extensions={"sgm", "sgml", "nrm", "ent","hdr", "hub", "dtd",
               "wml","vxml", "wml", "tld", "svg","xsl", "ecf", "jnlp", "xsd", "resx"} },
 { Lang="fsharp", Extensions={"fs","fsx"} },
 { Lang="informix", Extensions={"4gl"} },
 { Lang="blitzbasic", Extensions={"bb"} },
 { Lang="innosetup", Extensions={"iss"} },
 { Lang="lotus", Extensions={"ls"} },
 { Lang="ascend", Extensions={"a4c"} },
 { Lang="actionscript", Extensions={"as"} },
 { Lang="express", Extensions={"exp"} },
 { Lang="haxe", Extensions={"hx"} },
 { Lang="pyrex", Extensions={"pyx"} },

 { Lang="abap4", Extensions={"abp"} },
 { Lang="csharp", Extensions={"cs"} },
 { Lang="interlis", Extensions={"ili"} },
 { Lang="logtalk", Extensions={"lgt"} },
 { Lang="matlab", Extensions={"m"} },
 { Lang="nsis", Extensions={"nsi"} },
 { Lang="bison", Extensions={"y"} },
 { Lang="squirrel", Extensions={"nut"} },
 { Lang="luban", Extensions={"lbn"} },
 { Lang="maya", Extensions={"mel"} },
 { Lang="nemerle", Extensions={"n"} },
 { Lang="paradox", Extensions={"sc"} },
 { Lang="netrexx", Extensions={"nrx"} },
 { Lang="clearbasic", Extensions={"cb"} },
 { Lang="graphviz", Extensions={"dot"} },
 { Lang="small", Extensions={"sma"} },
 { Lang="autoit", Extensions={"au3"} },
 { Lang="chill", Extensions={"chl"} },
 { Lang="autohotkey", Extensions={"ahk"} },
 { Lang="fame", Extensions={"fame"} },
 { Lang="modelica", Extensions={"mo"} },
 { Lang="maple", Extensions={"mpl"} },
 { Lang="jasmin", Extensions={"j"} },
 { Lang="snobol", Extensions={"sno"} },
 { Lang="icon", Extensions={"icn"} },
 { Lang="felix", Extensions={"flx"} },
 { Lang="lindenscript", Extensions={"lsl"} },
 { Lang="lilypond", Extensions={"ly"} },
 { Lang="nasal", Extensions={"nas"} },
 { Lang="clean", Extensions={"icl"} },
 { Lang="assembler", Extensions={"asm"} },
 { Lang="bibtex", Extensions={"bib"} },
 { Lang="python", Extensions={"py"} },
 { Lang="txt", Extensions={"text"} },
 { Lang="n3", Extensions={"ttl", "nt"} },
 { Lang="biferno", Extensions={"bfr"} },
 { Lang="scilab", Extensions={"sci", "sce"} },

 { Lang="xml", Shebang=[[^\s*<\?xml\s+version=\"1\.0\"\s+[^(\?>)]*?>\s*$]] },
 { Lang="sh",  Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?([bd]ash|t?csh|[akz]?sh)]] },
 { Lang="make",Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?make]] },
 { Lang="awk", Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?[gnm]?awk]] },
 { Lang="perl",  Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?perl]] },
 { Lang="python",  Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?python]] },
 { Lang="ruby",  Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?ruby]] },
 { Lang="php", Shebang=[[^#!\s*(/usr)?(/local)?/bin/(env\s+)?php]] }
}
"""
import os
import sys
import unittest
import subprocess as SP
import mimetypes

from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.FileSystemXyz import BrowserWebView as Browser
import xyzpy.utilsXyz as UXYZ
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()


setTimer = True
deltaTime = 5000
verbose = False
verboseView = False
verboseEvent = False

testDir = os.path.realpath(
          os.path.split(
          os.path.realpath(__file__))[0])

_themesDir = "/usr/share/highlight/themes"

############################################
def isTextFileMimetypes(path):
  """if file acceptable for an ascii text editor"""
  mime = mimetypes.guess_type(path)
  if verbose: logger.info("mimetype: %s %s" % (mime,os.path.basename(path)))
  exclude = ["image"]
  mim0 = mime[0]
  if mim0 == None: return True #unknown mime as text editable (risky)
  for ex in exclude:
    if ex in mime[0]: return False
  return True

############################################
def isImageFileMimetypes(path):
  """is file image as tiff png jpeg etc..."""
  mime = mimetypes.guess_type(path)
  if verbose: logger.info("mimetype: %s %s" % (mime,os.path.basename(path)))
  mim0 = mime[0]
  if mim0 == None: return False #unknown mime as not image
  if "image" in mime[0]: return True
  return False

############################################
def isHtmlFileMimetypes(path):
  """is file html for direct display in web browser"""
  mime = mimetypes.guess_type(path)
  if verbose: logger.info("mimetype: %s %s" % (mime,os.path.basename(path)))
  mim0 = mime[0]
  if mim0 == None: return False #unknown mime as not html
  if "html" in mime[0]: return True
  return False

############################################
def isWebFileMimetypes(path):
  """is file acceptable for a web browser"""
  if isHtmlFileMimetypes(path): return True
  if isImageFileMimetypes(path): return True
  return False

############################################
class QdirViewTest(QtWidgets.QTreeView):

  def __init__(self, *args):
    super(QdirViewTest, self).__init__(*args)
    self._createActions()
    self._createContextMenus()

  def event(self, event):
    if verboseView: logger.info("QdirViewTest.event %s" % event)
    return super(QdirViewTest, self).event(event)

  def _createActions(self):
    self.expandActions = [
      self._createAction('Expand all', None, 'expand all items', self.expandAll) #, 'expand'),
    ]

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    if verbose: logger.info("QdirViewTest._createAction %s %s %s" % (self.__class__.__name__, Name, Icon))
    action = QtWidgets.QAction(Name, self)
    action.setIconVisibleInMenu(Icon != None)
    if Shortcut != None: 
      action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon != None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action

  def _createContextMenus(self):
    """each column have a different menu, virtual menus as example for inheritage"""
    menus = {}
    
    #outside column event
    menu = QtWidgets.QMenu("General", self)
    for action in self.expandActions: 
      menu.addAction(action)
    menus[-1] = menu
    
    #column 0 event
    #menu = QtWidgets.QMenu("Menu_0", self)
    #for action in #TODO: 
    #  menu.addAction(action)
    menus[0] = menu
    
    #column 1 event
    #menu = QtWidgets.QMenu("Menu_1", self)
    menus[1] = menu
    self.contextMenus = menus
    return

  def contextMenuEvent(self, event):
    """
    each column have a different menu, reals menus from user inheritage
    from user implementation of _createContextMenus
    """
    index = self.selectionModel().currentIndex()
    row = index.row()
    column = index.column()
    if verbose: logger.info("QdirViewTest.contextMenuEvent row %i column %i" % (row,column))
    
    pos = QtGui.QCursor.pos()
    theMenu = self.contextMenus[column]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    theMenu.exec_(pos)
    return

  def refreshFileView(self):
    selection = self.selectionModel()
    if selection.hasSelection():
      index = selection.currentIndex()
      self.clicked.emit(index) #will refresh fileView
      if verbose: logger.info("QdirViewTest.refreshFileView %s" % index)
    else:  # if there is no selection
      self.clicked.emit(None) #as root
      if verbose: logger.info("QdirViewTest.refreshFileView root")

############################################
class QfileViewTest(QtWidgets.QTreeView): 
  #QtWidgets.QTreeView or QtWidgets.QTableView
  
  helpSignal = QtCore.pyqtSignal()
  
  def __init__(self, *args):
    super(QfileViewTest, self).__init__(*args)
    #self.helpSignal.connect(self.on_help)
    #self.setToolTip("QfileViewTest")
    self._createActions()
    self._createContextMenus()
 
  def event(self, event):
    if verboseEvent: logger.info("QfileViewTest.event %s" % event)
    if type(event) == QtGui.QHelpEvent: 
      #print "QfileViewTest help event", event.type()
      self.helpSignal.emit()
    return super(QfileViewTest, self).event(event)

  def _createActions(self):
    editor = UXYZ.getEditor(basename=True)
    browser = UXYZ.getBrowser(basename=True)
    self.fileActions = [
      self._createAction('%s edit' % editor, None, 'edit by %s' % editor, self.editBy), #, 'expand'),
      self._createAction('%s view' % browser, None, 'view by %s' % browser, self.viewBy), #, 'expand'),
    ]

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    if verbose: logger.info("QdirViewTest._createAction %s %s %s" % (self.__class__.__name__, Name, Icon))
    action = QtWidgets.QAction(Name, self)
    action.setIconVisibleInMenu(Icon != None)
    if Shortcut != None: 
      action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon != None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action

  def _createContextMenus(self):
    """each column have a different menu, virtual menus as example for inheritage"""
    menus = {}
    
    #outside column event
    menu = QtWidgets.QMenu("General", self)
    for action in self.fileActions: 
      menu.addAction(action)
    menus[-1] = menu
    
    #column 0 event
    #menu = QtWidgets.QMenu("Menu_0", self)
    #for action in #TODO: 
    #  menu.addAction(action)
    menus[0] = menu
    
    #column 1 event
    #menu = QtWidgets.QMenu("Menu_1", self)
    menus[1] = menu
    self.contextMenus = menus
    return

  def contextMenuEvent(self, event):
    """
    each column have a different menu, reals menus from user inheritage
    from user implementation of _createContextMenus
    """
    index = self.selectionModel().currentIndex()
    row = index.row()
    column = index.column()
    self.currentPath = str(self.model().fileInfo(index).absoluteFilePath())
    #if verboseEvent: 
    logger.info("QfileViewTest.contextMenuEvent row %i column %i\n  '%s'" % (row,column,self.currentPath))
    
    pos = QtGui.QCursor.pos()
    theMenu = self.contextMenus[column]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    theMenu.exec_(pos)
    return

  def editBy(self):
    if verboseEvent: logger.info("editBy")
    if os.path.isdir(self.currentPath):
      QtWidgets.QMessageBox.warning(self,"warning","Problem is a directory:\n%s" % self.currentPath)
      return
    else:
      cmd = "%s %s &" % (UXYZ.getEditor(), self.currentPath)
      if verbose: logger.info("editBy stuff: '%s'" % cmd)
      process = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)

  def viewBy(self):
    if verboseEvent: logger.info("viewBy")
    cmd = "%s %s &" % (UXYZ.getBrowser(), self.currentPath)
    if verbose: logger.info("viewBy stuff: '%s'" % cmd)
    process = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)

  def xxon_help(self):
    index = self.selectionModel().currentIndex()
    row = self.selectionModel().currentIndex().row()
    column = self.selectionModel().currentIndex().column()
    model = self.model()
    #baseNode = self.selectedItems().getSelected[0]
    #for i in dir(self.selectionModel()): print "selectionModel",i
    #item = model.item(row,column)
    if verboseEvent: logger.info("QfileViewTest.on_help index %s %s %s" % (row,column,model))


############################################
class UserQFileSystemModel(QtWidgets.QFileSystemModel):

  def data(self, index, role):
    """to set tooltip"""
    if role == QtCore.Qt.ToolTipRole:
      #print "UserQFileSystemModel.data", index, role
      model = index.model() #PyQt5.QtWidgets.QFileSystemModel
      path = model.fileInfo(index).absoluteFilePath()
      return path #yet QtCore.QVariant(path)
    return super(UserQFileSystemModel, self).data(index, role)



  """do not work why??
  def setData(self, index, value, role):
    if True: #role == QtCore.Qt.ToolTipRole:
      print "UserQFileSystemModel.setData", index, role
    return super(UserQFileSystemModel, self).setData(index, value, role)"""


############################################
class QFileSystemModelDialog(QtWidgets.QWidget):
  """
  a widget dialog with a view for dirs and another one for current files
  as http://doc.qt.io/qt-5/model-view-programming.htmlh
  """
  def __init__(self, *args):
    super(QFileSystemModelDialog, self).__init__(*args)
    self.setWindowTitle("Files browser")
    layout = QtWidgets.QHBoxLayout()
    self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
    self.widget1 = QtWidgets.QWidget()
    self.widget2 = QtWidgets.QWidget()
    self.layout1 = QtWidgets.QVBoxLayout()
    self.layout2 = QtWidgets.QVBoxLayout()
    self.dirView = QdirViewTest()
    self.fileView = QfileViewTest()
    self.layout1.addWidget(self.dirView, 2)
    self.layout2.addWidget(self.fileView)
    self.widget1.setLayout(self.layout1)
    self.widget2.setLayout(self.layout2)
    self.splitter.addWidget(self.widget1)
    self.splitter.addWidget(self.widget2)
    layout.addWidget(self.splitter)
    self.setLayout(layout)

  def setDirRootPath(self, rootPath, filters=[]):
    """
    rootPath as "$PACKAGESPY_ROOT_DIR/toto"
    filters as ["*.py"]
    """
    realRootPath = self.getRealPath(rootPath)
    if not os.path.isdir(realRootPath):
      raise Exception("Inexisting directory: '%s'" % rootPath)
    self.originalRootPath = rootPath
    self.realRootPath = realRootPath 
    
    dirModel = UserQFileSystemModel() #QtWidgets.QFileSystemModel()
    fileModel = UserQFileSystemModel() #QtWidgets.QFileSystemModel()

    dirModel.setReadOnly(True)
    fileModel.setReadOnly(True)
    
    rootDir = dirModel.setRootPath(self.realRootPath)
    rootFile = fileModel.setRootPath(self.realRootPath)
    
    dirModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.AllDirs)
    fileModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.Files)
    fileModel.setNameFilters(filters)
    fileModel.setNameFilterDisables(False) #else only disable (grey)
    
    self.dirView.setModel(dirModel)
    self.fileView.setModel(fileModel)

    self.dirView.setRootIndex(rootDir)
    self.fileView.setRootIndex(rootFile)

    self.dirModel = dirModel
    self.fileModel = fileModel
    #print "fileModel:", id(self.fileModel)
    self.fileModel.userSelectedFiles = []

    self.dirView.hideColumn(1) #size
    self.dirView.hideColumn(2) #type
    self.dirView.hideColumn(3) #date modified
    
    self.fileView.header().setMinimumSectionSize(1)
    self.fileView.setColumnWidth(0, 200) #name
    self.fileView.setColumnWidth(1, 100) #size
    self.fileView.setColumnWidth(2, 1) #type
    self.fileView.setColumnWidth(3, 100) #date modified
    #self.fileView.hideColumn(2) #type
    
    self.fileView.setSortingEnabled(True)
    try:
      self.fileView.header().setSortIndicator(0, QtCore.Qt.AscendingOrder)
    except:
      #no header if self.fileView as QTableView
      pass

    #self.dirView.setRootIsDecorated(False) #useless
    #self.dirView.expandAll() #not working

    #self.setWindowTitle("Files browser (%s)" % self.realRootPath)
    self.setWindowTitle(self.realRootPath)

    self.theHeader=self.dirView.header()
    self.theHeader.setSectionsClickable(True)

    ##signals##
    self.theHeader.sectionClicked.connect(self.headerClikedEvent)
    self.theHeader.sectionDoubleClicked.connect(self.headerDoubleClickedEvent)

    #print "fileView.click",[i for i in dir(self.fileView) if "lick" in i]
    self.dirView.clicked.connect(self.on_dirView_clicked) #to refresh self.fileView
    self.fileView.clicked.connect(self.on_fileView_clicked) #to refresh self.webView/editView
    self.fileView.doubleClicked.connect(self.on_fileView_doubleClicked) #for example for futur
    #for i in dir(self.fileView): print self.fileView,i

  def on_dirView_clicked(self, index):
    if index == None: #root path
      path = self.realRootPath
    else:
      path = self.dirModel.fileInfo(index).absoluteFilePath()
    if verbose: logger.info("QFileSystemModelDialog.on_dirView_clicked slot on\n  '%s'" % path)
    self.setRootFile(path)

  def on_fileView_clicked(self, index):
    path = self.fileModel.fileInfo(index).absoluteFilePath()
    if verbose: logger.info("QFileSystemModelDialog.on_fileView_clicked slot on\n  '%s'" % path)
    if self.fileModel.isDir(index): #os.path.isdir(path): 
      self.setRootFile(path)
    else:
      #do nothing: edit/view if clicked see method inherited class UserQFileSystemModelDialog
      pass

  def on_fileView_doubleClicked(self, index):
    path = self.fileModel.fileInfo(index).absoluteFilePath()
    if verbose: logger.info("QFileSystemModelDialog.on_fileView_doubleClicked virtual slot on\n  '%s'" % path)

  def event(self, event):
    if verboseEvent: logger.info("QFileSystemModelDialog.event %s" % event)
    return super(QFileSystemModelDialog, self).event(event)

  def setRootFile(self, path):
    rootFile = self.fileModel.setRootPath(path)
    self.fileView.setRootIndex(rootFile)
    self.setWindowTitle(path)
    header = self.fileView.header()
    header.setStretchLastSection(False)
    self.fileView.resizeColumnToContents(0)
    #header.setSectionResizeMode(header.ResizeToContents)
    #header.setStretchLastSection(True)
    
  def headerClikedEvent(self, arg):
    if verboseEvent: logger.info("QFileSystemModelDialog.headerClikedEvent at column %i" % arg)
    self.setRootFile(self.realRootPath)
    self.dirView.selectionModel().clearSelection()

  def headerDoubleClickedEvent(self, arg):
    if verboseEvent: logger.info("QFileSystemModelDialog.headerDoubleClickedEvent at column %i" % arg)
    #self.setRootFile(self.realRootPath)
    #pos = QtGui.QCursor.pos()
    #theMenu = self.contextMenus[-1]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    #theMenu.exec_(pos)

  def getRealPath_obsolete(self, aPathWithEnvVar):
    """
    resolve env variable as $HOME/toto etc... 
    with subprocess shell interpretation of env var
    """
    cmd = "echo %s" % aPathWithEnvVar
    proc = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.rstrip()
    stderr = stderr.rstrip()
    #if verbose: print("stdout: '%s'\nstderr: '%s'" % (stdout, stderr))
    try:
      res = os.path.realpath(stdout)
    except:
      raise Exception("Problem shell interpretation of path: '%s'" % aPathWithEnvVar)
    
    if verbose: logger.info("getRealPath: '%s'->'%s'" % (aPathWithEnvVar, res))
    return res

  ################################################
  def getRealPath(self, aPathWithEnvVar):
    """
    resolve env variable as $HOME/toto etc...
    with expandvars
    """
    res = os.path.realpath(os.path.expandvars(aPathWithEnvVar))
    return res

  def toPathWithEnvVar(self, aPath):
    """
    returns a path with environ variable from original root path
    """
    aPathWithEnvVar = aPath.replace(self.realRootPath, self.originalRootPath)
    return aPathWithEnvVar


############################################
class UserItemDelegate(QtWidgets.QStyledItemDelegate):

  """example to change paint"""

  def paint(self, painter, option, index):

    modelItem = index.model() #PyQt5.QtWidgets.QFileSystemModel
    #print "delegate fileModel:", id(modelItem)
    #modelItem.parent()
    #treeItem = treeWidget.itemFromIndex(index)
    #print "treeItem",treeItem,index.row(),index().column()
    
    path = modelItem.fileInfo(index).absoluteFilePath()
    if path not in modelItem.userSelectedFiles:
      #print "not in userSelectedFiles",path
      textcolor = QtCore.Qt.black
    else:
      #print "in userSelectedFiles",path
      textcolor = QtCore.Qt.red
    
    painter.save()

    """# set background color
    painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))
    if option.state & QtWidgets.QStyle.State_Selected:
      painter.setBrush(QtGui.QBrush(QtCore.Qt.red))
    else:
      painter.setBrush(QtGui.QBrush(QtCore.Qt.white))
    painter.drawRect(option.rect)"""

    # set text color
    painter.setPen(QtGui.QPen(textcolor))
    value = index.data(QtCore.Qt.DisplayRole)
    painter.drawText(option.rect, QtCore.Qt.AlignLeft, value)
 
    painter.restore()
    return


############################################
class ComboBoxTheme(QtWidgets.QComboBox):
  def __init__(self, *args):
    super(ComboBoxTheme, self).__init__(*args)
    if os.path.isdir(_themesDir):
      res = os.listdir(_themesDir)
      res = [i.replace(".theme", "") for i in res if ".theme" in i]
      if verbose: logger.info("ComboBoxTheme themes:\n%s" % res)
      self.addItems(res)
      self.setValue(res[0])
    else:
      th = "no themes available"
      self.addItems([th])
      self.setValue(th)

  def setValue(self, value):
    if verbose: logger.info("ComboBoxTheme.setValue %s %s" % (str(value), type(value)))
    index = self.findText(str(value))
    self.setCurrentIndex(index)
    
  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.currentText())
  

############################################
class UserQStringListModel(QtCore.QStringListModel):
  #http://www.poketcode.com/en/pyqt4/item_views/index.html

  def data(self, index, role):
    """to set tooltip and display"""
    #print "UserQStringListModel.data", index, role

    if role == QtCore.Qt.ToolTipRole:
      #yet QtCore.QVariant(path)
      return super(UserQStringListModel, self).data(index, QtCore.Qt.EditRole)

    if role == QtCore.Qt.DisplayRole:
      pathname = super(UserQStringListModel, self).data(index, QtCore.Qt.EditRole)
      #print "pathname '%s'" % pathname
      path, name = os.path.split(pathname)
      return QtCore.QVariant(name)

    return super(UserQStringListModel, self).data(index, role)



############################################
class UserQListView(QtWidgets.QListView):

  """
  def contextMenuEvent(self, event): #TODO not used here
    index = self.selectionModel().currentIndex()
    row = index.row()
    column = index.column()
    print("UserListView.contextMenuEvent row %i column %i" % (row,column))
    
    pos = QtGui.QCursor.pos()
    menu = QtWidgets.QMenu("General", self)
    action = QtWidgets.QAction("remove", self)
    action.setIconVisibleInMenu(Icon != None)
    action.setToolTip("remove item in selected list")
    #action.triggered.connect(self.mainView.userSelectedFiles.remove(path))
    menu.addAction(action) 
    menu.exec_(pos)
    return
  """

  """
  def event(self, event):
    if verboseEvent: print "UserQListView.event",event.type(), event.MouseButtonDblClick
    #if event.type() == event.MouseButtonDblClick:
    #  print "!!!QEvent.MouseButtonDblClick"
    #  return
    #print "!!!QEvent.MouseButtonDblClick"
    return super(UserQListView, self).event(event)
   """

  """
  def mousePressEvent(self, event):
    if verboseEvent: print "UserQListView.mousePressEvent", event.type()
    return super(UserQListView, self).mousePressEvent(event)
  """

  def mouseDoubleClickEvent(self, event):
    """avoid user entry and modify view model"""
    if verboseEvent: logger.info("UserQListView.mouseDoubleClickEvent %s" % event.type())
    #return super(UserQListView, self).mouseDoubleClickEvent(event)
    return

############################################
class TextViewer(QtWidgets.QTextEdit): #or QtWidgets.QPlainTextEdit
  """
  viewer as read only syntax highliter
  """
  def __init__(self, *args):
    super(TextViewer, self).__init__(*args)
    self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.customContextMenuRequested.connect(self.contextMenuSlot)
    self.currentPath = None

  def _setNewContexMenu(self):
    self._menu = self.createStandardContextMenu()
    self._menu.addSeparator()
    quitAction = QtWidgets.QAction('%s edit' % UXYZ.getEditor(), self.editBy)
    quitAction.triggered.connect(self.test)
    self._menu.addAction(quitAction)

  def contextMenuSlot(self):
    self._setNewContexMenu()
    self._menu.exec_(QtGui.QCursor.pos())
  
  def editBy(self):
    if verboseEvent: logger.info("editBy")
    cmd = "%s %s &" % (UXYZ.getEditor(), self.currentPath)
    # print("begin %s '%s'" % (__file__, cmd))
    process = SP.Popen(cmd, shell=True)

############################################
class UserQFileSystemModelDialog(QFileSystemModelDialog):
  
  """
  QFileSystemModelDialog plus right third widget as textViewer & webViever of file
  """

  def __init__(self, *args):
    if verbose: logger.info("UserQFileSystemModelDialog.__init__")
    super(UserQFileSystemModelDialog, self).__init__(*args)
    self.splitter3 = QtWidgets.QSplitter(QtCore.Qt.Vertical)
    self.comboBoxTheme = ComboBoxTheme()
    self.editView = TextViewer() 
    # or QtWidgets.QTextEdit() #QtWidgets.QTextBrowser() #QtWidgets.QPlainTextEdit() #QtWidgets.QTextEdit()
    self.splitter3.addWidget(self.editView)
    self.webView = Browser() #QtWidgets.QTextBrowser()
    self.splitter3.addWidget(self.webView)
    self.splitter.addWidget(self.splitter3)
    self.layout1.addWidget(self.comboBoxTheme)
    self.selectedModel = UserQStringListModel([]) #QtCore.QStringListModel([])
    self.selectedView = UserQListView()
    self.selectedView.setModel(self.selectedModel)
    self.layout1.addWidget(self.selectedView, 1)

    self.currentEditViewFile = ""
    self.comboBoxTheme.currentIndexChanged.connect(self.on_changeTheme)
    self.fileView.setItemDelegate(UserItemDelegate(self))
    
    self.resize(1400,800)
    self.splitter.setSizes([200, 200, 1000])
    self.splitter3.setSizes([200, 0]) #hidden web view empty
    self.fileView.setColumnWidth(0, 200)

    self.editView.setReadOnly(True)
    self.editView.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse | QtCore.Qt.TextSelectableByKeyboard)


    #self.selectedView.doubleClicked.connect(self.on_selectedView_doubleClicked) #for example for futur
    
    #context menu
    self.selectedView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.selectedView.customContextMenuRequested.connect(self.openMenuSelectedView)

  def resetViewStretchForWeb(self):
    s = self.splitter3.sizes()
    if verbose: logger.info("reset sizes web",s)
    if s[1] < 20:
      self.splitter3.setSizes([100, 1000])

  def resetViewStretchForAscii(self):
    s = self.splitter3.sizes()
    if verbose: logger.info("reset sizes ascii",s)
    if s[0] < 20:
      self.splitter3.setSizes([1000, 100])

  def on_selectedView_doubleClicked(self):
    logger.info('on_selectedView_doubleClicked')

  def openMenuSelectedView(self, position):
    #https://wiki.python.org/moin/PyQt/Creating%20a%20context%20menu%20for%20a%20tree%20view
    if verbose: logger.info("openMenuSelectedView")
    index = self.selectedView.selectionModel().currentIndex()
    self.selectedModelIndexCurrent = index
    menu = QtWidgets.QMenu()
    action = QtWidgets.QAction("clear item", self)
    action.triggered.connect(self.clearOneSelectedModel)
    menu.addAction(action)
    action = QtWidgets.QAction("clear all", self)
    action.triggered.connect(self.clearAllSelectedModel)
    menu.addAction(action)
    menu.exec_(self.selectedView.viewport().mapToGlobal(position))

  def clearAllSelectedModel(self):
    """reset path in self.fileModel.selectedFiles"""
    if verboseEvent: logger.info("clearAllSelectedModel")
    self.fileModel.userSelectedFiles = []
    self.selectedModel.setStringList(self.fileModel.userSelectedFiles)
    self.dirView.refreshFileView()

  def clearOneSelectedModel(self):
    """reset path in self.fileModel.selectedFiles"""
    if verboseEvent: logger.info("clearOneSelectedModel %s" % self.selectedModelIndex.row())
    index = self.selectedModelIndexCurrent
    model = index.model()
    if model == None: return #model None in no selection
    path = model.data(index, QtCore.Qt.EditRole)
    self.fileModel.userSelectedFiles.remove(path)
    self.selectedModel.setStringList(self.fileModel.userSelectedFiles)
    self.dirView.refreshFileView()
    
  def on_changeTheme(self):
    if verboseEvent: logger.info("change theme %s" % self.comboBoxTheme.getValue())
    self.openFileHighlight(self.currentEditViewFile)
     
  def on_fileView_clicked(self, index):
    """edit the file clicked"""
    path = str(self.fileModel.fileInfo(index).absoluteFilePath())
    if os.path.isdir(path):
      self.setRootFile(path)
    else:
      if verboseEvent: logger.info("UserQFileSystemModelDialog.on_fileView_clicked view: '%s'" % path)
      #self.openFile(path)
      if isWebFileMimetypes(path):
        if verboseEvent: logger.info("Web %s" % mimetypes.guess_type(path))
        self.openBrowserFileHtml(path)
      if isTextFileMimetypes(path):
        if verboseEvent: logger.info("ascii %s" % mimetypes.guess_type(path))
        self.openFileHighlight(path)
      #self.openFileHtml(path)
      
  def on_fileView_doubleClicked(self, index):
    """set/reset path in self.fileModel.selectedFiles"""
    path = str(self.fileModel.fileInfo(index).absoluteFilePath())
    if path in self.fileModel.userSelectedFiles: #remove if present yet, else append
      self.fileModel.userSelectedFiles.remove(path)
    else:
      self.fileModel.userSelectedFiles.append(path)
    self.selectedModel.setStringList(self.fileModel.userSelectedFiles)  
    if verboseEvent: 
      logger.info("QFileSystemModelDialog.on_fileView_doubleClicked slot on:\n%s" % self.fileModel.userSelectedFiles)
    self.fileView.dataChanged(index,index) #refresh
    
  """
  def event(self, event):
    if verboseEvent: print "UserQFileSystemModelDialog.event",event
    if type(event) == QtGui.QHelpEvent: 
      if verboseEvent: print "help event", event.type()
    return super(QFileSystemModelDialog, self).event(event)
  """

  def openFile(self, nameFile):
    if nameFile == "": return #cancel
    """_, fileExtension = os.path.splitext(str(nameFile))
    if fileExtension in self.excludeFiles:"""
    if not isTextFileMimetypes(nameFile):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n%s" % nameFile)
      return
    if True:#try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        self.editView.clear()
        self.editView.setPlainText(open(nameFile, 'r').read())
        self.nameFile = nameFile
        self.currentEditViewFile = str(nameFile)
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    else: #except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)
    self.editView.currentFile = str(nameFile)

  def openFileHighlight(self, nameFile):
    if nameFile == "": return #cancel
    if not isTextFileMimetypes(str(nameFile)):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n%s" % nameFile)
      return
    if True: # try: TODO for debug
      if os.path.isfile(nameFile):
        #source accent is code python with #coding=utf-8
        self.editView.clear()
        cmd = "highlight -i %s --out-format=html --style %s --font-size=8" % (nameFile, self.comboBoxTheme.getValue())
        proc = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)
        res0, res1 = proc.communicate()
        res0 = UXYZ.toString(res0)
        res1 = UXYZ.toString(res1)
        #print "\n!!!!!!!!!!!!! stdout:\n",res0
        if res1 != "":
          logger.error("highlight:\n%s" % res1)
          self.openFile(nameFile)
        else:
          htmlStr = res0
          #self.editView.appendHtml(htmlStr)
          self.editView.insertHtml(htmlStr)
          self.nameFile = nameFile
          self.currentEditViewFile = str(nameFile)
          sb = self.editView.verticalScrollBar() #scroll to the top 
          sb.triggerAction(sb.SliderToMinimum)
          self.resetViewStretchForAscii()
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    else: # except: TODO for debug
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)
    self.currentEditViewFile = str(nameFile)
    self.editView.currentFile = str(nameFile)

  def openFileHtml(self, nameFile):
    if nameFile == "": return #cancel
    if not isTextFileMimetypes(str(nameFile)):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n%s" % nameFile)
      return
    if True:#try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        self.editView.clear()
        if True:
          self.editView.appendHtml(open(nameFile, 'r').read())
          self.nameFile = nameFile
          self.currentEditViewFile = str(nameFile)
          sb = self.editView.verticalScrollBar() #scroll to the top 
          sb.triggerAction(sb.SliderToMinimum)
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    else: #except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)
    self.currentEditViewFile = str(nameFile)
    self.editView.currentFile = str(nameFile)


  def openBrowserFileHtml(self, nameFile):
    #if not self.isWebview: return  #nothing
    if nameFile == "": return #cancel
    ok = isTextFileMimetypes(str(nameFile)) or isImageFileMimetypes(str(nameFile))
    if not ok:
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n%s" % nameFile)
      return
    if True:#try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        self.editView.clear()
        if True:
          self.webView.load(QtCore.QUrl(nameFile))
          self.nameFile = nameFile
          self.currentEditViewFile = str(nameFile)
          #sb = self.editView.verticalScrollBar() #scroll to the top 
          #sb.triggerAction(sb.SliderToMinimum)
          self.resetViewStretchForWeb()
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    else: #except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)
    self.currentEditViewFile = str(nameFile)



############################################
class TestCase(unittest.TestCase):
  
  def launchTimer(self, app, wid):
    if setTimer: 
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
    app.exec_()
  
  def xtest_100(self):
    app = OnceQApplication()
    currenDir = os.path.split(testDir)[0] #os.getcwd()
    treeView = QtWidgets.QTreeView()
    fileSystemModel = QtWidgets.QFileSystemModel(treeView)
    fileSystemModel.setReadOnly(True)
    root = fileSystemModel.setRootPath(currenDir)
    treeView.setModel(fileSystemModel)
    treeView.setRootIndex(root)
    treeView.setWindowTitle("test_100")
    treeView.show()
    self.launchTimer(app, treeView)

  def xtest_200(self):
    app = OnceQApplication()
    currenDir = os.path.split(testDir)[0] #os.getcwd()
    treeView = QfileViewTest()
    fileSystemModel = QtWidgets.QFileSystemModel(treeView)
    fileSystemModel.setReadOnly(True)
    root = fileSystemModel.setRootPath(currenDir)
    treeView.setModel(fileSystemModel)
    treeView.setRootIndex(root)
    treeView.setWindowTitle("test_200")
    treeView.show()
    self.launchTimer(app, treeView)

  def xtest_250(self):
    app = OnceQApplication()
    path = "$PACKAGESPY_ROOT_DIR"
    dialog = UserQFileSystemModelDialog()
    dialog.setDirRootPath(path)
    realPath = dialog.getRealPath(path)
    self.assertFalse("$" in realPath)
    newRealPath = os.path.join(realPath, "toto")
    newRealPathWithEnvVar = dialog.toPathWithEnvVar(newRealPath)
    self.assertTrue("$PACKAGESPY_ROOT_DIR" in newRealPathWithEnvVar)
    self.assertTrue("toto" in newRealPathWithEnvVar)

  def xtest_300(self):
    app = OnceQApplication()
    dialog = QFileSystemModelDialog()
    rootPath = "$PACKAGESPY_ROOT_DIR" 
    #os.path.realpath(os.path.join(testDir, '..', '..'))
    dialog.setDirRootPath(rootPath, filters=["*.py","*.c*","*.h*","*.C*"])
    dialog.show()
    self.launchTimer(app, dialog)

  def test_400(self):
    app = OnceQApplication()
    dialog = UserQFileSystemModelDialog()
    #path = "/volatile/wambeke/SAT4/prerequisites/uranie-3.6.0/FROM_nothing/share/uranie/html/index.html"
    #path = "helpWidgetUserQFileSystemModelDialog.html" TODO?
    #dialog.webView.load(QtCore.QUrl(path))
    #rootPath = os.path.join("$PACKAGESPY_ROOT_DIR", 'pythonAppliMatix')
    #rootPath = os.path.join("$URANIESYS", "share", "uranie")
    rootPath = os.path.join(testDir, "..", "..")
    #dialog.setDirRootPath(rootPath, filters=["*.py","*.c*","*.h*","*.C*"])
    dialog.setDirRootPath(rootPath, filters=[])
    f = dialog.fileView
    #for i in dir(f): print "dialog.fileView",i
    dialog.show()
    self.launchTimer(app, dialog)

if __name__ == '__main__':
  verbose = False
  setTimer = False #user call and wait
  unittest.main()
  pass
