

### PACKAGESPY tests

Launch tests with bash script `.../scripts/AllTestLauncher.sh`.
This script launch `.../pythonAppliMatix/unittestpy/AllTestLauncher.py`.

#### Prerequisites

- As *standalone mode*, these python tests need Python3 PyQt5 etc... environment,
  *for example* use [miniconda](https://docs.conda.io/en/latest/miniconda.html)
  ```bash
  conda activate yourPythonEnvironment   # a python3 PyQt5 numpy matplotlib etc.
  ```

- ... or ... as *SALOME context mode*
  ```bash
  cd .../SALOME-9.9.0       # example where there is SALOME
  ./salome context
  ```

- ... or ... as *MATIX context mode*
  ```bash
  cd .../MATIX_V30          # example where there is MATIX
  ./matix context
  ```

#### Launch tests

- Use AllTestLauncher.sh
  ```bash
  cd ${PACKAGESPY_ROOT_DIR}/pythonAppliMatix  # ... or ... cd .../PACKAGESPY/pythonAppliMatix
  # all tests
  ./scripts/AllTestLauncherPACKAGESPY.sh . "test_???_*.py"
  # one test
  ./scripts/AllTestLauncherPACKAGESPY.sh . "test_174_*.py"
  ```

- Or directly use of AllTestLauncher.py
  ```bash
  cd ${PACKAGESPY_ROOT_DIR}/pythonAppliMatix
  # all tests
  ./unittestpy/AllTestLauncher.py -r pythonAppliMatix -p "test_???_*.py"
  # one test
  ./unittestpy/AllTestLauncher.py -r pythonAppliMatix -p "test_174_*.py"
  ```


### AllTestLauncher.py

See usage:

```bash
${PACKAGESPY_ROOT_DIR}/pythonAppliMatix/unittestpy/AllTestLauncher.py --help

    usage: AllTestLauncher.py [-h] [-d] [-v int] [-r dirPath] [-p filePattern]
                              [-t outputType] [-n dirName]

    launch All python tests

    optional arguments:
      -h, --help            show this help message and exit
      -d, --debug           set debug mode, more verbose
      -v int, --verbosity int
                            set verbosity of unittests [0|1|2...]
      -r dirPath, --rootPath dirPath
                            dir name with absolute or relative path stand for root
                            directory of recursive searching unittest python files
      -p filePattern, --pattern filePattern
                            file pattern for unittest files
                            ['test_*.py'|'*Test.py'...]
      -t outputType, --type outputType
                            type of output: ['std'(standart ascii)|'xml'|'html']
      -n dirName, --name dirName
                            (only for type xml) name of directory output:
                            ['test_reports'|...]. If name = 'stdout' then all-in-
                            one xml output at 'sys.stdout'. For pipe redirection:
                            '>> AllTestLauncher.py -t xml -n stdout > tmp.xml'
```
