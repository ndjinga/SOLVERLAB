#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

r'''
Easy interface using inspect and pydoc

WARNING: getfullargspec(method) #only python 3
'''

#see https://tomassetti.me/python-reflection-how-to-list-modules-and-inspect-functions/

import os
import pprint as PP
import inspect
import pydoc
from collections import namedtuple

verbose = False

########################################################
class Inspect(object):
  def __init__(self, anObject):
    self._object = anObject
    self._class = self._object.__class__.__name__
    self._delimiter = "\n############ %s:\n%s\n############"
    self._DefaultArgSpec = namedtuple('DefaultArgSpec', 'has_default default_value')


  def getWithDelimiter(self, res):
    """for comprehensive successive prints"""
    return self._delimiter % (self._class, res)

  def get_inspect_getargspec(self, method):
    return PP.pformat(inspect.getargspec(method))
    try:
      res = inspect.getfullargspec(method) #python 3
      return res
    except:
      res = inspect.getargspec(method)

  def _get_default_arg(self, args, defaults, arg_index):
    """
    Method that determines if an argument has default value or not,
    and if yes what is the default value for the argument

    :param args: array of arguments, eg: ['first_arg', 'second_arg', 'third_arg']
    :param defaults: array of default values, eg: (42, 'something')
    :param arg_index: index of the argument in the argument array for which,
    this function checks if a default value exists or not. And if default value
    exists it would return the default value. Example argument: 1
    :return: Tuple of whether there is a default or not, and if yes the default
    value, eg: for index 2 i.e. for "second_arg" this function returns (True, 42)
    """
    if not defaults:
        return self._DefaultArgSpec(False, None)

    args_with_no_defaults = len(args) - len(defaults)

    if arg_index < args_with_no_defaults:
        return self._DefaultArgSpec(False, None)
    else:
        value = defaults[arg_index - args_with_no_defaults]
        if (type(value) is str):
            value = '"%s"' % value
        return self._DefaultArgSpec(True, value)

  def get_inspect_signature(self, method):
    """
    Given a function, it returns a string that pretty much looks how the
    function signature would be written in python.

    :param method: a python method
    :return: A string similar describing the pythong method signature.
    eg: "my_method(first_argArg, second_arg=42, third_arg='something')"
    """

    # The return value of ArgSpec is a bit weird, as the list of arguments and
    # list of defaults are returned in separate array.
    # eg: ArgSpec(args=['first_arg', 'second_arg', 'third_arg'],
    # varargs=None, keywords=None, defaults=(42, 'something'))
    argspec = inspect.getargspec(method)
    arg_index=0
    args = []

    # Use the args and defaults array returned by argspec and find out
    # which arguments has default
    for arg in argspec.args:
        default_arg = self._get_default_arg(argspec.args, argspec.defaults, arg_index)
        if default_arg.has_default:
            args.append("%s=%s" % (arg, default_arg.default_value))
        else:
            args.append(arg)
        arg_index += 1
    return "%s(%s)" % (method.__name__, ", ".join(args))

  def get_pydoc_render_doc(self):
    """
    removed the boldface sequences with pydoc.plain
    python -c 'import pydoc; print pydoc.render_doc(pydoc)'
    """
    return pydoc.plain(pydoc.render_doc(self._object, title="Documentation %s"))
 
  def get_prefixedMethodsDoc(self, prefix=""):
    """
    example usage for ControllerIpc:
    b = ControllerIpc()
    print(Inspect(b).get_PrefixedMethodsDoc("Ipc"))
    """
    # inspect.getdoc reindent docstrings
    lg = len(prefix)
    if lg > 0:
      Methods = [ (name, self.get_inspect_signature(method), inspect.getdoc(method)) \
         for name, method in inspect.getmembers(self._object, predicate=inspect.ismethod) \
         if prefix == name[:lg] ]
    else:
      Methods = [ (name, self.get_inspect_signature(method), inspect.getdoc(method)) \
         for name, method in inspect.getmembers(self._object, predicate=inspect.ismethod) ]
    
    #print( PP.pformat(Methods) )

    res = ""
    for name, signature, docstr in Methods:
      if docstr == None: docstr = "No documentation"
      indentedDoc = "\n".join(["    " + line for line in docstr.split("\n")])
      res += "\n  METHOD: %s\n%s\n" % (signature, indentedDoc)
    return res


