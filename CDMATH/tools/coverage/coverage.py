import os
import subprocess
import sys
from optparse import OptionParser


def launchCoverage(cdmath_build, src_cdmath, step, module_cdmath, coverage):
    directory_object = os.path.join(cdmath_build, module_cdmath, "src",
                                    "CMakeFiles", module_cdmath + ".dir")
    fi1 = os.path.join(src_cdmath, 'pre_requis/MEDCoupling')
    fi2 = os.path.join(src_cdmath, 'pre_requis/MEDLoader')
    if coverage == "python":
        extention = "py"
        pass
    else:
        extention = "cpp"
        pass

    if not os.path.isdir(directory_object):
        sys.stderr.write('Error:  cdmath directory objects ' +
                         directory_object + ' does not exist \n')
        sys.exit(1)
        pass

    if step == "a":
        proc = subprocess.Popen("lcov --zerocounters --directory " +
                                directory_object, shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
        proc = subprocess.Popen("lcov --capture --initial --directory " + directory_object + " --output-file=" +
                                module_cdmath + "_init_" + extention + ".info", shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
        pass
    else:
        proc = subprocess.Popen("lcov --capture --directory " + directory_object + " --output-file=" +
                                module_cdmath + "_tests_" + extention + ".info", shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
        proc = subprocess.Popen("lcov -a " + module_cdmath + "_init_" + extention + ".info -a " + module_cdmath + "_tests_" + extention +
                                ".info -o " + module_cdmath + "_tests_all_" + extention + ".info", shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
        proc = subprocess.Popen("lcov --remove " + module_cdmath + "_tests_all_" + extention + ".info '/usr/include/c++/4.7/x86_64-linux-gnu/bits/*' '/usr/include/*' '/usr/include/c++/4.7/*' '" +
                                fi1 + "/*' '" + fi2 + "/*' --output " + module_cdmath + "_" + extention + ".info", shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
        proc = subprocess.Popen('genhtml --title \"CDMATH COVERAGE TESTS\" --show-details --highlight --legend --output-directory html_' +
                                extention + ' ' + module_cdmath + '_' + extention + '.info', shell=True, stdout=sys.stdout, stderr=sys.stderr)
        so, se = proc.communicate()
    return


def main(argv):
    cdmath_build, src_cdmath, step, module_cdmath, coverage = options_script(argv)
    launchCoverage(cdmath_build, src_cdmath, step, module_cdmath, coverage)
    return


def options_script(argv):

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-b", "--build-dir", dest="build_cdmath",
                      metavar="<build_cdmath>",
                      help="choose your build directory.")

    parser.add_option("-d", "--src-dir", dest="src_cdmath",
                      metavar="<src_cdmath>",
                      help="choose your src directory.")

    parser.add_option("-m", "--module", dest="module_cdmath",
                      metavar="<module_cdmath>",
                      help="choose your cdmath module.")

    parser.add_option("-s", "--step", dest="step", choices=["b", "a"],
                      metavar="<step>",
                      help="choose your step. ")

    parser.add_option("-c", "--coverage", dest="coverage", choices=["python", "cpp"],
                      metavar="<coverage>",
                      help="choose your coverage type (C++ or python). ")

#    parser.set_defaults(library_type="STATIC")

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)
        pass

    build_cdmath = os.path.expanduser(options.build_cdmath)
    build_cdmath = os.path.expandvars(build_cdmath)
    build_cdmath = os.path.abspath(build_cdmath)
    if not os.path.isdir(build_cdmath):
        sys.stderr.write('Error: build cdmath directory ' + build_cdmath + ' does not exist \n')
        sys.exit(1)
        pass

    src_cdmath = os.path.expanduser(options.src_cdmath)
    src_cdmath = os.path.expandvars(src_cdmath)
    src_cdmath = os.path.abspath(src_cdmath)
    if not os.path.isdir(src_cdmath):
        sys.stderr.write('Error: src cdmath directory ' + src_cdmath + ' does not exist \n')
        sys.exit(1)
        pass

    module_cdmath = options.module_cdmath
    coverage = options.coverage
    step = options.step
    return build_cdmath, src_cdmath, step, module_cdmath, coverage

if __name__ == "__main__":
    main(sys.argv[1:])
