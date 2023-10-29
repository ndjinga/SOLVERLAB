#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Logger for packagespy using logging package

used by solverlabGui and other iradinaGUI as possibly outside SALOME Desktop.

| Define one logger with one handler on stdout and one handler on file for production.
| Define another one logger for unittest.
|
| see: http://sametmax.com/ecrire-des-logs-en-python
|
| Define two LoggerSvl instances in solverlabGui, no more need.
|   - _loggerDefault as production/development logger
|   - _loggerUnittest as unittest logger
|
| see use of handlers of _loggerDefault for
| log console and log files xml, txt
|
| console handler:
|   - info() : no format
|   - error() warning() trace() debug() etc. :
|      formatted indented on multi lines messages using handlers
|
| file handlers:
|   - info() error() warning() trace() debug() etc. :
|      formatted indented on multi lines messages using handlers
|
| WARNING:
|  log step and log trace are present on stdout console or log file
|  following level handlers settings
"""

import os
import sys
import time
import logging as LOGI
from logging.handlers import BufferingHandler
import pprint as PP

_verbose = False # True
_name = "loggingSvl"
_loggerDefaultName = 'SvlDefaultLogger'
_loggerUnittestName = 'SvlUnittestLogger'

_loggerCurrent = []

_MISSED_GETENV = []

levels = {
  "CRITICAL": LOGI.CRITICAL,
  "ERROR": LOGI.ERROR,
  "WARNING": LOGI.WARNING,
  "INFO": LOGI.INFO,
  #"STEP": _STEP, # step level is just below INFO
  #"TRACE": _TRACE, # trace level is just below STEP,
  "DEBUG": LOGI.DEBUG
}

_levelNames = {
  LOGI.CRITICAL: "CRITICAL",
  LOGI.ERROR: "ERROR",
  LOGI.WARNING: "WARNING",
  LOGI.INFO: "INFO",
  #LOGI.STEP: "STEP",
  #LOGI.TRACE: "TRACE",
  LOGI.DEBUG: "DEBUG"
}


if _verbose:
  print("CRITICAL  %s" % LOGI.CRITICAL)
  print("WARNING  %s" % LOGI.WARN)
  print("INFO  %s" % LOGI.INFO)

_knownLevels = "CRITICAL ERROR WARNING INFO DEBUG".split()
_knownLevelsStr = "[%s]" % "|".join(_knownLevels)

#################################################################
# utilities methods
#################################################################

def filterLevel(aLevel):
  """
  filter levels logging values from firsts characters levels.
  No case sensitive.

  | example:
  | 'i' -> 'INFO'
  | 'cRiT' -> 'CRITICAL'
  """
  aLev = aLevel.upper()
  knownLevels = _knownLevels
  maxLen = max([len(i) for i in knownLevels])
  for i in range(maxLen):
    for lev in knownLevels:
      if aLev == lev[:i]:
        # DBG.write("filterLevel", "%s -> %s" % (aLevel, lev))
        return lev
  msg = "Unknown level '%s', accepted are: %s" % (aLev, ", ".join(knownLevels))
  raise Exception(msg)
  # or signal problem but message debug
  # DBG.write("Problem filterLevel", msg, True)
  # return "CRITICAL"

def indent(msg, nb, car=" "):
  """indent nb car (spaces) multi lines message except first one"""
  s = msg.split("\n")
  res = ("\n"+car*nb).join(s)
  return res

def indentUnittest(msg, prefix=" | "):
  """
  indent multi lines message except first one with prefix.
  prefix default is designed for less spaces for size logs files
  and keep logs human eye readable
  """
  s = msg.split("\n")
  res = ("\n" + prefix).join(s)
  return res

def log(msg, force=False):
  """elementary log when no logging.Logger yet"""
  prefix = "---- %s.log: " % _name
  nb = len(prefix)
  if _verbose or force:
    print(prefix + indent(msg, nb))

# just for debug info where is import logging
log("import logging on %s" % LOGI.__file__)


def getStrDirLogger(logger):
  """
  Returns multi line string for logger description, with dir(logger).
  Used for debug
  """
  lgr = logger # shortcut
  msg = "%s(name=%s, dateLogger=%s):\n%s\n"
  cName = lgr.__class__.__name__
  res = msg % (cName, lgr.name, lgr.dateLogger, PP.pformat(dir(lgr)))
  return res

def getStrHandler(handler):
  """
  Returns one line string for handler description
  (as inexisting __repr__)
  to avoid create inherited classe(s) handler
  """
  h = handler # shortcut
  msg = "%s(name=%s)"
  cName = h.__class__.__name__
  res = msg % (cName, h.get_name())
  return res

def getStrShort(msg):
  """Returns short string for msg (as first caracters without line feed"""
  # log("getStrShort " + str(msg), True)
  res = msg.replace("\n", "//")[0:30]
  return res

def getStrLogRecord(logRecord):
  """
  Returns one line string for simple logging LogRecord description
  """
  import solverlabpy.coloringSvl as COLS # needs sys.path.insert(0, svldir)
  msg = "LogRecord(level='%s', msg='%s...')"
  shortMsg = getStrShort(logRecord.msg)
  levelName = COLS.cleanColors(logRecord.levelname).replace(" ", "")
  res = msg % (levelName, shortMsg)
  return res

def getListOfStrLogRecord(listOfLogRecord):
  """
  Returns one line string for logging LogRecord description
  """
  res = [getStrLogRecord(l) for l in listOfLogRecord]
  return res

#################################################################
# solverlabGui logger classes
#################################################################

try:
    str
    _unicode = True
except NameError:
    _unicode = False

def getMessage(self):
    """
    modified from logging.__init__.LogRecord.getMessage,
    better message on format error
    Return the message for this LogRecord.

    Return the message for this LogRecord after merging any user-supplied
    arguments with the message.
    """
    if not _unicode: #if no unicode support...
        msg = str(self.msg)
    else:
        msg = self.msg
        if not isinstance(msg, str):
            try:
                msg = str(self.msg)
            except UnicodeError:
                msg = self.msg      #Defer encoding till later
    if self.args:
        try: # better message on format error
          msg = msg % self.args
        except Exception as e:
          msg = "ERROR: %s with args %s" % (msg, PP.pformat(self.args))
          log(msg, True)
    return msg

LOGI.LogRecord.getMessage = getMessage # better message if error

#################################################################
class LoggerSvl(LOGI.Logger):
  """
  Inherited class logging.Logger for logger solverlabGui

  see: /usr/lib64/python2.7/logging/__init__.py etc.
  """

  def __init__(self, name, level=LOGI.INFO):
    """
    Initialize the logger with a name and an optional level.
    """
    super(LoggerSvl, self).__init__(name, level)
    self.dateLogger = "NoDateLogger"
    self.dateHour = None # datehour of main command
    self.isClosed = False
    self.idCommandHandlers = 0 # incremented, 0 for main command 1, 2, etc. for micro command

    self._logStep_header = ["No log step header ..."] # one line message for steps
    self._loggerParamiko = None # set only if solverlab jobs, useless for now
    self._fileParamiko = None # set when main command set solverlab logger

  def logStep_begin(self, header, step=""):
    """
    initialize for main handler (tty as stdout)
    a one line message for steps (...of compilation for example)
    as no return line message with logger.info() level

    | example:
    | 'header ... first temporary step message ...'
    | 'header ... etc ...'         (on same line)
    | 'header ... OK'              (on same line)
    """
    self._logStep_header.append(header)
    self.logStep(step)

  def logStep(self, step):
    """
    current logger.info() step as
    'header ... etc ...'
    """
    header = self._logStep_header[-1]
    if type(step) == str:
      self.info("<RC>%s %s ..." % (header, step))
      return
    elif step.isOk(): # as ReturnCode type step
      self.info("<RC>%s <OK> ..." % header)
    else:
      self.info("<RC>%s %s <KO> ..." % (header, step.getWhy()))

  def logStep_end(self, step, tab=None):
    """
    last logger.info() step as
    'header ... OK' or 'header ... KO'
    """
    import solverlabpy.utilsSvl as UTS

    header = self._logStep_header[-1]
    if tab is None:
      if type(step) == str:
        self.info("<RC>%s %s" % (header, step))
      elif step.isOk(): # as ReturnCode type
        self.info("<RC>%s <OK>" % header)
      else:
        self.info("<RC>%s <%s>" % (header, step.getStatus()))
      # pop as end
      if len(self._logStep_header) > 1:
        self._logStep_header.pop()
      else:
        self.error("Something wrong for logStep")
      return

    else:
      if type(step) == str:
        stepTab = UTS.tabColor(tab, header, 0, step)
        self.info("<RC>%s" % stepTab)
      elif step.isOk(): # as ReturnCode type
        stepTab = UTS.tabColor(tab, header, 0, "<OK>")
        self.info("<RC>%s" % stepTab)
      else:
        stepTab = UTS.tabColor(tab, header, 0, "<%s>: %s" % (step.getStatus(), step.getWhy()))
        self.info("<RC>%s" % stepTab)
      # pop as end
      if len(self._logStep_header) > 1:
        self._logStep_header.pop()
      else:
        self.error("Something wrong for logStep")
      return

  def getMainCommandHandler(self):
    """
    returns handler for colored stdout console/terminal
    for human user eye solverlabGUI outputs
    """
    for h in self.handlers:
      if h.idCommandHandlers == 0:
        return h
    return None

  def close(self):
    """
    final stuff for logger, done at end solverlabGui
    flushed and closed xml files have to be not overriden/appended
    """
    if self.isClosed:
      raise Exception("logger closed yet: %s" % self)
    log("close stuff logger %s" % self) # getStrDirLogger(self)
    for handl in list(self.handlers):  # get original list
      log("close stuff handler %s" % getStrHandler(handl))
      handl.close() # Tidy up any resources used by the handler.
      self.removeHandler(handl)
    # todo etc
    self.isClosed = True # done at end solverlabGUI, flushed closed xml files.
    return

  def __repr__(self):
    """one line string representation"""
    msg = "%s(name=%s, dateLogger=%s, handlers=%s)"
    cName = self.__class__.__name__
    h = [getStrHandler(h) for h in self.handlers]
    h = "[" + ", ".join(h) + "]"
    res = msg % (cName, self.name, self.dateLogger, h)
    return res

  def xx_isEnabledFor(self, level):
    """
    Is this logger enabled for level 'level'?
    currently not modified from logging.Logger class,
    here only for call log debug.
    """
    log("logger %s isEnabledFor %i>=%i" % (self.name, level, self.getEffectiveLevel()))
    if self.manager.disable >= level:
        return 0
    return level >= self.getEffectiveLevel()

  def setFileHandlerForCommand(self, cmdParent, cmdInstance):
    """
    add file handler to logger to set log files
    for a solverlabGui command.
    when command is known from pyconf/config instance

    | Example:
    | log files names for command prepare
    | with micro commands clean/source/patch
    |   ~/LOGS/20180510_140606_prepare_lenovo.xml
    |   ~/LOGS/OUT/20180510_140606_prepare_lenovo.txt
    |   ~/LOGS/micro_20180510_140607_clean_lenovo.xml
    |   ~/LOGS/OUT/micro_20180510_140607_clean_lenovo.txt
    |   etc.
    """
    import solverlabpy.utilsSvl as UTS # needs sys.path.insert(0, svldir)

    logger = self
    config = cmdInstance.getConfig()

    #import debogpy.debug as DBG # avoid cross import
    log("setFileHandler %s" % logger)
    log("setFileHandler config\n%s" % PP.pformat(dict(config.VARS)))
    log("setFileHandler TODO set log_dir config.LOCAL.log_dir")

    log_dir = config.LOCAL.log_dir # files xml
    log_dir_out = os.path.join(log_dir, "OUT") # files txt
    log_dir_jobs = os.path.join(log_dir, "JOBS") # files txt
    UTS.ensure_path_exists(log_dir)
    UTS.ensure_path_exists(log_dir_out)
    if self.idCommandHandlers == 0:
      datehour = config.VARS.datehour
      self.dateHour = datehour # save dateHour for micro commands
    else:
      datehour = self.dateHour # micro commands have same datehour for naming files

    cmd = config.VARS.command
    fullNameCmd = cmdInstance.getFullNameStr()
    hostname = config.VARS.hostname
    nameFileXml = "%s_%03i_%s_%s.xml" % (datehour, self.idCommandHandlers, cmd, hostname)
    nameFileTxt = "%s_%03i_%s_%s.txt" % (datehour, self.idCommandHandlers, cmd, hostname)
    fileXml = os.path.join(log_dir, nameFileXml)
    fileTxt = os.path.join(log_dir_out, nameFileTxt)

    # precaution
    lastCmd = cmdInstance.getFullNameList()[-1]
    if cmd != lastCmd:
      msg = "setFileHandler '%s' command name incoherency in config '%s'" % (fullNameCmd, cmd)
      logger.critical(msg)

    nbhandl = len(logger.handlers) # number of active current handlers

    if self.idCommandHandlers == 0: # first main command
      log("setFileHandler '%s' main command (id=%i)" % (fullNameCmd, self.idCommandHandlers))

      ################################
      # Logging vers file xml
      handler = XmlHandler(3000) # no many log outputs in memory
      handler.setLevel(LOGI.STEP)
      handler.set_name(nameFileXml)
      handler.set_target_file(fileXml)
      handler.set_config(config)
      handler.idCommandHandlers = self.idCommandHandlers

      fmt = '%(asctime)s :: %(levelname)-8s :: %(message)s'
      formatter = FileXmlFormatter(fmt, "%y-%m-%d %H:%M:%S")

      handler.setFormatter(formatter)
      logger.addHandler(handler)

      ################################
      # Logging vers file txt
      handler = LOGI.FileHandler(fileTxt)
      handler.setLevel(LOGI.TRACE)
      handler.set_name(nameFileTxt)
      handler.idCommandHandlers = self.idCommandHandlers

      fmt = '%(asctime)s :: %(levelname)-8s :: %(message)s'
      formatter = FileTxtFormatter(fmt, "%y-%m-%d %H:%M:%S")

      handler.setFormatter(formatter)
      logger.addHandler(handler)


    elif self.idCommandHandlers > 0: # secondary micro command
      log("TODO setFileHandler '%s' micro command (id=%i)" % (fullNameCmd, self.idCommandHandlers))

      ################################
      # Logging vers file xml
      handler = XmlHandler(3000) # no many log outputs in memory
      handler.setLevel(LOGI.STEP)
      handler.set_name(nameFileXml)
      handler.set_target_file(fileXml)
      handler.set_config(config)
      handler.idCommandHandlers = self.idCommandHandlers

      fmt = '%(asctime)s :: %(levelname)-8s :: %(message)s'
      formatter = FileXmlFormatter(fmt, "%y-%m-%d %H:%M:%S")

      handler.setFormatter(formatter)
      logger.addHandler(handler)

      ################################
      # Logging vers file txt
      handler = LOGI.FileHandler(fileTxt)
      handler.setLevel(LOGI.TRACE)
      handler.set_name(nameFileTxt)
      handler.idCommandHandlers = self.idCommandHandlers

      fmt = '%(asctime)s :: %(levelname)-8s :: %(message)s'
      formatter = FileTxtFormatter(fmt, "%y-%m-%d %H:%M:%S")

      handler.setFormatter(formatter)
      logger.addHandler(handler)

    cmdInstance.setIdCommandHandlers(self.idCommandHandlers)
    newLink = self.initLinkForCommand(cmdParent, cmdInstance)
    newLink.setAuthAttr("cmd_name", cmd)
    newLink.setAuthAttr("log_file_name", nameFileXml)

    self.idCommandHandlers += 1
    log("setFileHandler %s" % logger)
    return self.idCommandHandlers

  def setLevelMainHandler(self, level):
    for handl in list(self.handlers): # get main handler
      if handl.idCommandHandlers == 0:
        log("setLevelMainHandler %s" % level)
        handl.setLevel(level)
        return
    raise Exception("main handler not found for level %s" % level)

  def closeFileHandlerForCommand(self, cmdInstance):
    for handl in list(self.handlers): # get original list
      try: # may be foreign handlers without idCommandHandlers attribute
        if handl.idCommandHandlers == cmdInstance._idCommandHandlers:
          nbini = len(self.handlers)
          log("close stuff handler %s" % getStrHandler(handl))
          handl.close() # Tidy up any resources used by the handler.
          self.removeHandler(handl)
          if len(self.handlers) != nbini-1:
            self.critical("Unexpected len(logger.handlers)=%i" %  len(self.handlers))
      except:
        self.warning("Existing logger handler without idCommandHandlers attribute %s" % str(handl))
        pass

  def initLinkForCommand(self, cmdParent, cmdNew):
    import solverlabpy.linksXml as LIXML
    newLink = LIXML.appendLinkForCommand(cmdParent, cmdNew)
    return newLink

  def testNoReturn(self):
    """test when message ending '...' and level info then no return mode"""
    testNoReturn(self) # module method

  def getenv(self, name, default=""):
    """
    signal only one time missing env var
    returns empty string in this case, not None
    """
    res = os.getenv(name, default=None)
    if res is None:
      if name not in _MISSED_GETENV:  # warning problem log only one time
        _MISSED_GETENV.append(name)
        self.warning("missed environment variable '%s', default '%s'" % (name, default))
        res = default  # avoid raise error, problems in os.path.join() AttributeError: 'NoneType' object has no attribute 'endswith'
    # print("log getenv '%s'" % res)
    return res


#################################################################
class DefaultFormatter(LOGI.Formatter):

  # to set color prefix, problem with indent format as
  _ColorLevelname = {
    "DEBUG": "<green>",
    "TRACE": "<green>",
    "STEP": "<green>",
    "INFO":  "<green>",
    "WARNING": "<red>",
    "ERROR": "<yellow>",
    "CRITICAL": "<yellow>",
  }

  def format(self, record):
    import solverlabpy.coloringSvl as COLS # needs sys.path.insert(0, svldir)
    if _verbose:
      import debogpy.debug as DBG # avoid cross import
      DBG.write("DefaultFormatter.format", "%s: %s..." % (record.levelname, record.msg[0:20]), True)
    record.levelname = self.setColorLevelname(record.levelname)
    # INFO messages could be NOT prefixed
    if "xxxINFO" in record.levelname: # but xxxINFO as every times prefixed
      res = str(record.msg)
    else:
      res = indent(super(DefaultFormatter, self).format(record), 12)
    res = COLS.toColor(res)
    # print "begin_mes'%s'end_mess" % res
    return res

  def setColorLevelname(self, levelname):
    """
    set color implies color special characters and
    tabulate levelname length of string
    """
    color = self._ColorLevelname[levelname]
    res = color + levelname + "<reset>"
    nb = len(levelname)
    res = res + " "*(8-nb) # 8 as len("CRITICAL")
    # log("setColorLevelname'%s'" % res)
    return res


#################################################################
class UnittestFormatter(LOGI.Formatter):
  def format(self, record):
    import solverlabpy.coloringSvl as COLS # needs sys.path.insert(0, svldir)
    # print "", record.levelname #type(record), dir(record)
    # nb = len("2018-03-17 12:15:41 :: INFO     :: ")
    res = super(UnittestFormatter, self).format(record)
    res = indentUnittest(res)
    return COLS.toColor(res)


#################################################################
class FileTxtFormatter(LOGI.Formatter):
  def format(self, record):
    import solverlabpy.coloringSvl as COLS # needs sys.path.insert(0, svldir)
    # print "", record.levelname #type(record), dir(record)
    # nb = len("2018-03-17 12:15:41 :: INFO     :: ")
    res = super(FileTxtFormatter, self).format(record)
    res = indentUnittest(res)
    return COLS.cleanColors(res)


#################################################################
class FileXmlFormatter(LOGI.Formatter):
  def format(self, record):
    import solverlabpy.coloringSvl as COLS # needs sys.path.insert(0, svldir)
    # print "", record.levelname #type(record), dir(record)
    # nb = len("2018-03-17 12:15:41 :: INFO     :: ")
    res = super(FileXmlFormatter, self).format(record)
    res = indentUnittest(res)
    return COLS.cleanColors(res)


#################################################################
class UnittestStream(object):
  """
  write my stream class
  only write and flush are used for the streaming

  | https://docs.python.org/2/library/logging.handlers.html
  | https://stackoverflow.com/questions/31999627/storing-logger-messages-in-a-string
  """
  def __init__(self):
    self._logs = ''

  def getLogs(self):
    return self._logs

  def getLogsAndClear(self):
    res = self._logs
    self._logs = ''
    return res

  def write(self, astr):
    """final method called when message is logged"""
    # log("UnittestStream.write('%s')" % astr) # for debug ...
    self._logs += astr

  def flush(self):
    pass

  def __str__(self):
    return self._logs

#################################################################
class StreamHandlerSvl(LOGI.StreamHandler):
    """
    A handler class which writes logging records, appropriately formatted,
    to a stream. Note that this class does not close the stream, as
    sys.stdout or sys.stderr may be used.

    from logging.StreamHandler class,
    modified for 'no return' mode line if '...' at end of record message
    """
    def isLastRecordHaveNoReturn(self):
        """
        to memorize if last info record is 'no return' mode (as ending '...')
        avoid define inherited __init__
        """
        if not hasattr(self, "lastRecordHaveNoReturn"):
          self.lastRecordHaveNoReturn = False
        return self.lastRecordHaveNoReturn

    def isNeedFirstReturn(self, record):
        """
        'no return' mode  valid only if 2 consecutives info messages
        if not, needs insert return line BEFORE (warning, debug, or other)
        current record message
        """
        if not self.isLastRecordHaveNoReturn():
          return False
        if record.levelno == LOGI.INFO:
          return False # case is 2 consecutives info messages (OK for continuity)
        return True # case need insert return line BEFORE message (KO for continuity)

    def emit(self, record):
        """
        Emit a record.

        If a formatter is specified, it is used to format the record.
        The record is then written to the stream with a trailing newline.  If
        exception information is present, it is formatted using
        traceback.print_exception and appended to the stream.  If the stream
        has an 'encoding' attribute, it is used to determine how to do the
        output to the stream.
        """
        try:
            msg = self.format(record)
            stream = self.stream
            if msg[-3:] == "..." and record.levelno == LOGI.INFO:
              fs = '%s'
              ufs = '%s'
              self.lastRecordHaveNoReturn = True
            else:
              if self.isNeedFirstReturn(record):
                fs = '\n%s\n'
                ufs = '\n%s\n'
                self.lastRecordHaveNoReturn = False
              else:
                fs = '%s\n'
                ufs = '%s\n'
            if not _unicode: #if no unicode support...
                stream.write(fs % msg)
            else:
                try:
                    if (isinstance(msg, str) and
                        getattr(stream, 'encoding', None)):
                        # ufs = u'%s\n'
                        try:
                            stream.write(ufs % msg)
                        except UnicodeEncodeError:
                            #Printing to terminals sometimes fails. For example,
                            #with an encoding of 'cp1251', the above write will
                            #work if written to a stream opened or wrapped by
                            #the codecs module, but fail when writing to a
                            #terminal even when the codepage is set to cp1251.
                            #An extra encoding step seems to be needed.
                            stream.write((ufs % msg).encode(stream.encoding))
                    else:
                        stream.write(fs % msg)
                except UnicodeError:
                    stream.write(fs % msg.encode("UTF-8"))
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

#################################################################
class XmlHandler(BufferingHandler):
  """
  log outputs in memory as BufferingHandler.
  Write ElementTree in file and flush are done once
  when method close is called, to generate xml file.

  | atts = {
  |   "fileName": xml file name of micro command
  |   "command": cmd, # 'compile' or 'prepare' etc.
  |   "passed": res, # 'O' or '1'
  |   "launchedCommand" : fullcmd, # 'compile TOTO -etc'
  |   }
  |
  | see: https://docs.python.org/2/library/logging.handlers.html
  """
  def __init__(self, capacity):
    super(XmlHandler, self).__init__(capacity)
    self._target_file = None
    self._config = None
    self._log_field = "Uninitialized log"
    self._final_fields = {} # node attributes
    self.isClosed = False # precaution as write file done yet
    self.idCommandHandlers = None # have to be set later to know links

  def set_target_file(self, filename):
    """
    filename is file name xml with path
    supposedly non existing, no overwrite accepted
    """
    if os.path.exists(filename):
      msg = "XmlHandler target file %s existing yet" % filename
      raise Exception(msg)
    self._target_file = filename

  def set_config(self, config):
    """
    config is supposedly non existing, no overwrite accepted
    """
    if self._config is not None:
      msg = "XmlHandler target config existing yet"
      raise Exception(msg)
    self._config = config

  def close(self):
    """
    prepare ElementTree from existing logs and write xml file

    warning: avoid solverlabGUI logging message in logger close phase
    """
    targetFile = self._target_file
    config = self._config

    # log("dir(XmlHandler)\n" + PP.pformat(dir(self)), True)

    if self.isClosed:
      msg = "XmlHandler target file %s closed yet" % targetFile
      log(msg, True) #avoid solverlabGUI logging message in logger close phase
      return # avoid overwrite

    if os.path.exists(targetFile):
      msg = "XmlHandler target file %s existing yet" % targetFile
      log(msg, True) #avoid solverlabGUI logging message in logger close phase
      return # avoid overwrite
    """
    else: # for debug
      msg = "XmlHandler target file %s NOT existing yet" % targetFile
      log(msg, True) #avoid solverlabGUI logging message in logger close phase
    """

    # TODO for debug
    log("XmlHandler to xml file\n%s" % PP.pformat(getListOfStrLogRecord(self.buffer)))

    self._log_field = self.createLogField()

    xmlFile = XMLMGR.XmlLogFile(targetFile, "solverlabGUIcommand")
    xmlFile.put_initial_fields(config)
    xmlFile.put_log_field(self._log_field)
    xmlFile.put_links_fields(self.idCommandHandlers)
    xmlFile.put_final_fields(self._final_fields)
    xmlFile.write_tree(stylesheet = "command.xsl") # xml complete closed file
    xmlFile.dump_config(config) # create pyconf file in the log directory

    self.isClosed = True # precaution to not override xml closed file
    # zaps the buffer to empty as parent class
    super(XmlHandler, self).close() # n.b. extract handler from logger

  def createLogFieldFromScrath(self):
    """
    prepare formatted string from self.buffer LogRecord for xml 'Log' node
    local format
    """
    import solverlabpy.coloringSvl as COLS
    res = ""
    for lr in self.buffer:
      fmt = "%s :: %s\n"
      if lr.levelno != LOGI.INFO:
        levelName = COLS.cleanColors(lr.levelname).replace(" ", "")
        msg = COLS.cleanColors(lr.msg)
        res += fmt % (levelName, msg)
    if res == "":
      res = "Empty log"
    return res

  def createLogField(self):
    """
    prepare formatted string from self.buffer LogRecord for xml 'Log' node
    using handler formatter
    """
    import solverlabpy.coloringSvl as COLS
    fmtr = self.formatter
    res = ""
    for lr in self.buffer:
      if lr.levelno != LOGI.INFO: #skip info level and supposed no debug present
        res += fmtr.format(lr) + "\n"
    if res == "":
      res = "Empty log"
    return COLS.cleanColors(res)


#################################################################
# methods to define two LoggerSvl instances in solverlabGui,
# no more need
#################################################################
def initLoggerAsDefault(logger, fmt=None, level=None):
  """
  init logger as prefixed message and indented message if multi line
  exept info() outed 'as it' without any format.
  level could be modified during execution
  """
  log("initLoggerAsDefault name=%s\nfmt='%s' level='%s'" % (logger.name, fmt, level))
  handler = StreamHandlerSvl(sys.stdout) # Logging vers console
  handler.set_name(logger.name + "_console")
  if fmt is not None:
    # formatter = LOGI.Formatter(fmt, "%Y-%m-%d %H:%M:%S")
    formatter = DefaultFormatter(fmt, "%y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
  handler.idCommandHandlers = 0
  logger.addHandler(handler)
  # as RootLogger is level WARNING
  # my logger is not notset but low, handlers needs setlevel greater
  logger.setLevel(LOGI.DEBUG)
  # import src/debug as DBG
  # tmp = (logger.getEffectiveLevel(), LOGI.NOTSET, logger.level, logger.parent.level)
  # DBG.write("logger levels tmp, True)
  if level is not None: # level could be modified during execution
    handler.setLevel(level) # on screen log as user wants
  else:
    handler.setLevel(LOGI.INFO) # on screen no log step, which are in xml files
  return


def initLoggerAsUnittest(logger, fmt=None, level=None):
  """
  init logger as silent on stdout/stderr
  used for retrieve messages in memory for post execution unittest
  https://docs.python.org/2/library/logging.handlers.html
  """
  log("initLoggerAsUnittest name=%s\nfmt='%s' level='%s'" % (logger.name, fmt, level))
  stream = UnittestStream()
  handler = LOGI.StreamHandler(stream) # Logging vers stream
  handler.set_name(logger.name + "_unittest")
  if fmt is not None:
    # formatter = LOGI.Formatter(fmt, "%Y-%m-%d %H:%M:%S")
    formatter = UnittestFormatter(fmt, "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
  handler.idCommandHandlers = 0
  logger.addHandler(handler)
  logger.stream = stream
  logger.getLogs = stream.getLogs
  logger.getLogsAndClear = stream.getLogsAndClear
  if level is not None:
    logger.setLevel(level)
  else:
    logger.setLevel(LOGI.DEBUG)


def _getDefaultLogger():
  """internal use only"""
  log("getDefaultLogger %s" % _loggerDefaultName)
  # case multithread may be problem as not LOGI._acquireLock()
  previousClass = LOGI._loggerClass
  LOGI.setLoggerClass(LoggerSvl) # to get LoggerSvl instance with trace etc.
  res = LOGI.getLogger(_loggerDefaultName)
  LOGI.setLoggerClass(previousClass)
  return res


def _getUnittestLogger():
  """internal use only"""
  log("getUnittestLogger %s" % _loggerUnittestName)
  # case multithread may be problem as not LOGI._acquireLock()
  previousClass = LOGI._loggerClass
  LOGI.setLoggerClass(LoggerSvl) # to get LoggerSvl instance with trace etc.
  res = LOGI.getLogger(_loggerUnittestName)
  LOGI.setLoggerClass(previousClass)
  return res


def getDefaultLogger():
  """official method to get the only one instance of Default Logger"""
  res = _getDefaultLogger()
  _loggerCurrent.append(res)
  return res


def getLogger():
  return getDefaultLogger()


def getUnittestLogger():
  """official method to get the only one instance of Unittest Logger"""
  res = _getUnittestLogger()
  _loggerCurrent.append(res)
  return res


def getCurrentLogger():
  """
  get one of _loggerDefault or _loggerUnittest as first created by
  getDefautLogger() or getUnittestLogger():
  """
  nb = len(_loggerCurrent)
  if nb != 1:
    log("WARNING: no _loggerCurrent for now, create DefaultLogger", force=True)
    res = getDefaultLogger()
    return res
  return _loggerCurrent[-1]


#################################################################
def getenv(name, default=""):
  return getLogger().getenv(name, default)


#################################################################
# small tests as demonstration, see unittest also
#################################################################
def testLogger_1(logger):
  """small test"""
  # print getStrDirLogger(logger)
  logger.debug('test logger debug')
  logger.info('test logger info')
  logger.warning('test logger warning')
  logger.error('test logger error')
  logger.critical('test logger critical')
  logger.info('\ntest logger info: no indent\n- second line\n- third line\n')
  logger.warning('test logger warning:\n- second line\n- third line')


def testNoReturn(logger):
  """test when message ending '...' and level info then no return mode"""
  logger.warning('BEGIN test NoReturn 1')
  logger.info('test no return here 0...')
  logger.info('have to continue here 1...')
  logger.info('have to continue here 2')
  logger.info('END test NoReturn 1')
  logger.warning('BEGIN test NoReturn 2')
  logger.info('test no return here 0...')
  logger.warning('have NOT to continue here 1...')
  logger.info('have NOT to continue here 2')
  logger.info('END test NoReturn 2')

def pushLevel(level):
  #obsolete
  return

def popLevel():
  #obsolete
  return None

def testMain():
  print("\n**** DEFAULT logger")
  logdef = getDefaultLogger()
  # use of setColorLevelname <color>...<reset>, so do not use %(levelname)-8s
  initLoggerAsDefault(logdef, '%(levelname)s :: %(message)s', level=LOGI.DEBUG)
  testLogger_1(logdef)
  print("\n**** DEFAULT logger NoReturn")
  logdef.testNoReturn()

  print("\n**** UNITTEST logger")
  loguni = getUnittestLogger()
  initLoggerAsUnittest(loguni, '%(asctime)s :: %(levelname)-8s :: %(message)s', level=LOGI.DEBUG)
  testLogger_1(loguni) # is silent
  # log("loguni.getLogs():\n%s" % loguni.getLogs())
  print("loguni.streamUnittest:\n%s" % loguni.getLogs())
  print("\n**** UNITTEST logger NoReturn")
  loguni.testNoReturn()
  print("loguni.streamUnittest:\n%s" % loguni.getLogs())

  from colorama import Fore as FG
  from colorama import Style as ST
  print("this is unconditionally %scolored in green%s !!!" % (FG.GREEN, ST.RESET_ALL))

  import solverlabpy.utilsSvl as UTS
  import solverlabpy.coloringSvl as COLS
  print("\n1234567890123456789012345678901234567890123456789012345678901234567890")
  print(UTS.tabColor(10, "0", 10, "1", 10, "2", 10, "3", 10, "4", 10, "5", 10, "6"))
  print(UTS.tabColor(20, "1 tabulated", 15, "21 OK here", 10, "36 OK end"))
  print(COLS.toColor(UTS.tabColor(20, "1 <green>tabulated<reset>", 15, "21 <OK> <info>here<reset>", 10, "36 <OK> end")))
  print('toColor(20, "1 <green>tabulated<reset>", 15, "21 <OK> <info>here<reset>", 10, "36 <OK> end")')

  print("\n1234567890123456789012345678901234567890123456789012345678901234567890")
  print(UTS.tabColor(10, "0", 10, "1", 10, "2", 10, "3", 10, "4", 10, "5", 10, "6"))
  print(UTS.tabColor(-20, "tabulated 20", -15, "OK here 35", -10, "OK end 45"))
  print(COLS.toColor(UTS.tabColor(-20, "<green>tabulated<reset> 20", -15, "<OK> <info>here<reset> 35", -10, "<OK> end 45")))
  print('toColor(-20, "<green>tabulated<reset> 20", -15, "<OK> <info>here<reset> 35", -10, "<OK> end 45")')



#################################################################
# in production, or not (if __main__)
#################################################################
if __name__ == "__main__":
  # for example, not in production
  # get path to solverlabGui sources
  svldir = os.path.dirname(os.path.dirname(__file__))
  # Make the src & commands package accessible from all code
  sys.path.insert(0, svldir)
  testMain()
  # here we have sys.exit()
else:
  # in production
  # get two LoggerSvl instance used in solverlabGui, no more needed.
  _loggerDefault = _getDefaultLogger()
  _loggerUnittest = _getUnittestLogger()
  # initLoggerAsDefault(_loggerDefault, '%(levelname)s :: %(name)s :: %(message)s')
  initLoggerAsDefault(_loggerDefault, '%(levelname)s : %(message)s')
  initLoggerAsUnittest(_loggerUnittest, '%(asctime)s :: %(levelname)s :: %(message)s')
