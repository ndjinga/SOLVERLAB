Build on an online build service
================================

Open Build Service
------------------
Here is how to create `cdmath_0.3-1.dsc` and `cdmath_0.3-1.tar.gz`, necessary for the Open Build Service (provided by SUSE):
 * be on a Debian-based system,
 * copy the source directory of CDMATH (`cdmath_src` probably) and name this copy `cdmath-0.3`,
 * run `dpkg-source -b cdmath-0.3/`.

Launchpad.net
-------------
Here is how to create files, necessary for the Launchpad.net build service (provided by Canonical):
 * be on a Debian-based system,
 * copy the source directory of CDMATH (`cdmath_src` probably) and name this copy `cdmath_0.3`,
 * compress `cdmath_0.3` to `cdmath_0.3-1.orig.tar.gz`,
 * move to the `cdmath_0.3` directory,
 * run `debuild -S`,
 * upload with `dput ppa:cdmath/cdmath cdmath_0.3-1trusty_source.changes`.
