%module cdmath

%pythoncode %{
def FieldIadd(self,*args):
    import _cdmath
    return _cdmath.Field____iadd___(self,self, *args)
def FieldIsub(self,*args):
    import _cdmath
    return _cdmath.Field____isub___(self,self, *args)
def FieldImul(self,*args):
    import _cdmath
    return _cdmath.Field____imul___(self,self, *args)
def FieldIdiv(self,*args):
    import _cdmath
    return _cdmath.Field____idiv___(self,self, *args)

def MatrixIadd(self,*args):
    import _cdmath
    return _cdmath.Matrix____iadd___(self,self, *args)
def MatrixIsub(self,*args):
    import _cdmath
    return _cdmath.Matrix____isub___(self,self, *args)
def MatrixImul(self,*args):
    import _cdmath
    return _cdmath.Matrix____imul___(self,self, *args)
def MatrixIdiv(self,*args):
    import _cdmath
    return _cdmath.Matrix____idiv___(self,self, *args)

def VectorIadd(self,*args):
    import _cdmath
    return _cdmath.Vector____iadd___(self,self, *args)
def VectorIsub(self,*args):
    import _cdmath
    return _cdmath.Vector____isub___(self,self, *args)
def VectorImul(self,*args):
    import _cdmath
    return _cdmath.Vector____imul___(self,self, *args)
def VectorIdiv(self,*args):
    import _cdmath
    return _cdmath.Vector____idiv___(self,self, *args)
%}

%{
#include "Matrix.hxx"
#include "Vector.hxx"
#include "Point.hxx"
#include "Node.hxx"
#include "Cell.hxx"
#include "Face.hxx"
#include "Mesh.hxx"
#include "Field.hxx"

#include <sstream>

%}

%feature("autodoc", "1");
%feature("docstring");

#define MEDCOUPLING_EXPORT

%include "MEDCouplingCommon.i"
%include std_string.i

%include "GenericMatrix.hxx"
%include "Matrix.hxx"
%include "Vector.hxx"
%include "Point.hxx"
%include "Node.hxx"
%include "Cell.hxx"
%include "Face.hxx"
%include "Mesh.hxx"
%include "Field.hxx"

%extend Vector
{
  Vector __sub__(const Vector& f)
  {
    return (*self)-f;
  }

  Vector __add__(const Vector& f)
  {
    return (*self)+f;
  }

  double __mul__(const Vector& f)
  {
    return (*self)*f;
  }

  Vector __mul__(double value)
  {
    return (*self)*value;
  }

  Vector __mul__(int value)
  {
    return (*self)*value;
  }

  Vector __div__(double value)
  {
    return (*self)/value;
  }

  Vector __div__(int value)
  {
    return (*self)/value;
  }

  double __getitem__(int i)
  {
    const double& tmp=(*self)(i);
    return double(tmp);
  }

  void __setitem__(int i, double val)
  {
    (*self)(i)=val;
  }

  PyObject *___iadd___(PyObject *trueSelf, const Vector& f)
  {
    (*self)+=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___iadd___(PyObject *trueSelf, double value)
  {
    (*self)+=value;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___isub___(PyObject *trueSelf, const Vector& f)
  {
    (*self)-=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___isub___(PyObject *trueSelf, double value)
  {
    (*self)-=value;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___imul___(PyObject *trueSelf, double val)
  {
    (*self)*=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___idiv___(PyObject *trueSelf, double val)
  {
    (*self)/=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  std::string __repr__()
  {
     std::ostringstream oss;
     oss << "Number of Values : " << (*self).getNumberOfRows() << std::endl;
     oss << "Data content : " << std::endl;
     for (int i=0;i<(*self).getNumberOfRows();i++)
         oss << "# " << i << " : " << (*self)(i) << std::endl;
     return oss.str();
  }
}

%extend Matrix
{
  Matrix __sub__(const Matrix& f)
  {
    return (*self)-f;
  }

  Matrix __add__(const Matrix& f)
  {
    return (*self)+f;
  }

  Matrix __mul__(const Matrix& f)
  {
    return (*self)*f;
  }

  Matrix __mul__(double value)
  {
    return (*self)*value;
  }

  Matrix __mul__(int value)
  {
    return (*self)*value;
  }

  Matrix __div__(double value)
  {
    return (*self)/value;
  }

  Matrix __div__(int value)
  {
    return (*self)/value;
  }

  double __getitem__(PyObject* args)
  {
      int sz(PyTuple_Size(args));
      if(sz!=2)
         throw std::out_of_range("number of arguments!=2");
      PyObject *pyI(PyTuple_GetItem(args,0)),*pyJ(PyTuple_GetItem(args,1));
      int i(PyInt_AS_LONG(pyI)),j(PyInt_AS_LONG(pyJ));
      const double& tmp=(*self)(i,j);
      return double(tmp);
  }

  void __setitem__(PyObject* args, double val)
  {
    int sz(PyTuple_Size(args));
    if(sz!=2)
      throw std::out_of_range("number of arguments!=3");
    PyObject *pyI(PyTuple_GetItem(args,0)),*pyJ(PyTuple_GetItem(args,1));
    int i(PyInt_AS_LONG(pyI)),j(PyInt_AS_LONG(pyJ));
    (*self)(i,j)=val;
  }

  PyObject *___iadd___(PyObject *trueSelf, const Matrix& f)
  {
    (*self)+=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___isub___(PyObject *trueSelf, const Matrix& f)
  {
    (*self)-=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___imul___(PyObject *trueSelf, double val)
  {
    (*self)*=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___imul___(PyObject *trueSelf, const Matrix& f)
  {
    (*self)*=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___idiv___(PyObject *trueSelf, double val)
  {
    (*self)/=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  std::string __repr__()
  {
     std::ostringstream oss;
     oss << "Number of Values : " << (*self).getNumberOfRows() << std::endl;
     oss << "Data content : " << std::endl;
     for (int i=0;i<(*self).getNumberOfRows();i++)
         for (int j=0;j<(*self).getNumberOfColumns();j++)
         oss << "# (" << i << "," << j << ") : " << (*self)(i,j) << std::endl;
     return oss.str();
  }
}

%extend Point
{
  double __getitem__(int i)
  {
    const double& tmp=(*self)[i];
    return double(tmp);
  }

  void __setitem__(int i, double val)
  {
    (*self)[i]=val;
  }
}

%extend Mesh
{
  std::string __repr__()
  {
     std::ostringstream oss;
     oss << self->getMEDCouplingMesh()->simpleRepr();
     return oss.str();
    }

}

%extend Field
{
  Field __sub__(const Field& f)
  {
    return (*self)-f;
  }

  Field __add__(const Field& f)
  {
    return (*self)+f;
  }

  Field __mul__(double value)
  {
    return (*self)*value;
  }

  Field __mul__(int value)
  {
    return (*self)*value;
  }

  Field __div__(double value)
  {
    return (*self)/value;
  }

  Field __div__(int value)
  {
    return (*self)/value;
  }

  PyObject *___iadd___(PyObject *trueSelf, const Field& f)
  {
    (*self)+=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___iadd___(PyObject *trueSelf, double value)
  {
    (*self)+=value;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___isub___(PyObject *trueSelf, const Field& f)
  {
    (*self)-=f;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___isub___(PyObject *trueSelf, double value)
  {
    (*self)-=value;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___imul___(PyObject *trueSelf, double val)
  {
    (*self)*=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *___idiv___(PyObject *trueSelf, double val)
  {
    (*self)/=val;
    Py_XINCREF(trueSelf);
    return trueSelf;
  }

  PyObject *getValues()
  {
     const double *vals=self->getValues();
     return convertDblArrToPyList(vals,self->getNumberOfElements()*self->getNumberOfComponents());
  }

  double __getitem__(PyObject *value)
  {
    int sz=PyTuple_Size(value);
    PyObject *ob1=PyTuple_GetItem(value,0);
    PyObject *ob2=PyTuple_GetItem(value,1);
    int i=(int)PyInt_AS_LONG(ob1);
    int j=(int)PyInt_AS_LONG(ob2);
    const double& tmp=(*self)(i,j);
    return double(tmp);
  }

  double __getitem__(int i)
  {
    const double& tmp=(*self)(i);
    return double(tmp);
  }

  void __setitem__(int i, double val)
  {
    (*self)(i)=val;
  }

  void __setitem__(PyObject *ob, double value)
  {
    int sz=PyTuple_Size(ob);
    PyObject *ob1=PyTuple_GetItem(ob,0);
    PyObject *ob2=PyTuple_GetItem(ob,1);
    int i=(int)PyInt_AS_LONG(ob1);
    int j=(int)PyInt_AS_LONG(ob2);
    (*self)(i,j)=value;
  }

  std::string __repr__()
  {
     std::ostringstream oss;
     oss << self->getField()->simpleRepr();
     return oss.str();
    }
}

%pythoncode %{
Field.__iadd__=FieldIadd
Field.__isub__=FieldIsub
Field.__imul__=FieldImul
Field.__idiv__=FieldIdiv

Matrix.__iadd__=MatrixIadd
Matrix.__isub__=MatrixIsub
Matrix.__imul__=MatrixImul
Matrix.__idiv__=MatrixIdiv

Vector.__iadd__=VectorIadd
Vector.__isub__=VectorIsub
Vector.__imul__=VectorImul
Vector.__idiv__=VectorIdiv
%}

namespace std {
 %template(VectorVectorDouble) vector< vector<double> >;
};
