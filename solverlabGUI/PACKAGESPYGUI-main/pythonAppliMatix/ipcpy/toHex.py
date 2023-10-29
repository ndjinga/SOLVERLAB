#!/usr/bin/env python
# coding: utf-8 

"""
simple hexadecimal viewer utility
"""

#####################################################
def toHexAll(aStr):
  res = " ".join(i.encode("hex") for i in aStr)
  return res

#####################################################
def toHex(aStr):
  """return line every 10 hexa representations"""
  res = toHexAll(aStr).split(" ")
  rres = "\n"
  icou = 0
  while icou < len(res):
    hexList = res[icou: icou+10]
    rres += "{0:>3}: {1:<30}// {2:<12}\n".format(str(icou), " ".join(hexList), toStr(hexList))
    icou += 10
  return rres

#####################################################
def toStr(hexList, default="^"):
  """filter special characters as default '^'"""
  #chr(int(hex_data[x:y], 16)
  res = ""
  for i in hexList:
    try:
      ii = int(i, 16)
      if ii >= 32 and ii <= 126:
        res += chr(ii)
      else:
        res += default #"."
    except:
      res += default #"."
  return "'" + res + "'"
  #rres = str(bytearray.fromhex(hexStr).decode())
  return "'" + rres + "'"
  
#####################################################
def test_toHex():
  aStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!!!!hexa xFE:\xFEoctal007:\007tab:\t     CR:\r      LF:\n      NULL:\x00"
  print("INFO: test_toHex for:", type(aStr)) #, "\n" , toHexAll(aStr)
  print("INFO: test_toHex result:",toHex(aStr))

if __name__ == "__main__":
  test_toHex()

