import glob
# import pdb
import string
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import dpath
import VTK_routines
import PV_routines
import MED_routines
# import pprint as pp
import logging
from matplotlib.pyplot import xlabel
import __main__


# TODO : to be replaced with a proper log system
try:
    LOG_FILE_NAME = os.path.splitext(__main__.__file__)[0] + '.log'
except:
    LOG_FILE_NAME = 'resultsBaseProcess.log'
logging.basicConfig(filename=LOG_FILE_NAME, level=logging.DEBUG)
# logging.basicConfig(filename=__main__.__file__+'resDb.log', level=logging.WARNING)

# TODO : shouldn't we use join??


def create_file_list(workDir, ext=None, fname_glob=None):
    expanded_path = os.path.expanduser(workDir)
    abs_path = os.path.abspath(expanded_path)
    # TODO : we should test if directory exists
    if fname_glob is None:
        pattern_for_files = abs_path + '/*' + ext
    elif ext is None:
        pattern_for_files = abs_path + '/' + fname_glob
    file_list = glob.glob(pattern_for_files)
    return file_list

# TODO : we should use regexp instead?



# TODO : the new tag system should be implemented
class Result:
    def __init__(self, file_name,
                 variable,
                 dimension,
                 instant=None,
                 step=None,
                 solver=None,
                 tags=[],
                 code=None,
                 is_stationnary=True
                 ):
        self.file_name = os.path.abspath(file_name)
        absolute_root_name, ext = os.path.splitext(self.file_name)
        root_name = os.path.basename(absolute_root_name)
        self.root_name = root_name
        self.file_extension = ext
        self.tags = set()
        self.tags.update(set(tags))
        self.code = code
        self.solver = solver
        self.variable = variable
        self.instant = instant
        self.step = step
        self.is_stationnary = is_stationnary
        self.dim = dimension
        self.default_probe = self.create_default_probe()
#        self.checkResFieldIsValid() # useless
        for k in ['file_extension', 'code', 'variable',
                  'instant', 'step', 'is_stationnary', 'dim', 'default_probe', 'solver']:
            self.tags.add('required:{}={}'.format(k, self.__dict__[k]))

    def __repr__(self):
        txt = ""
        for k, v in self.__dict__.items():
            txt += "\t{}={}\n".format(k, v)
        return txt

# does this really make sense?
    @property
    # TODO : can be replaced with a prop
    def setTags(self, newTags):
        self.tags + newTags

class Probe:
    def __init__(self, file_name,variable_name,
                 name,
                 line=None,
                 ExtractData=None,
                 tags=[]
                 ):
        self.file_name = file_name
        self.variable_name=variable_name
        self.probe_name = name
        self.line = line
        self.extract_data = ExtractData
        self._tags = set(['required:probe_id={}'.format(self)])
        self._tags.update(tags)

    @property
    def tags(self):
        return self._tags.union(self.field_result.tags)

    @property
    def master_tag(self):
        return string.join(self._tags.union(self.field_result.tags), sep=";")

    @property
    def data(self):
        return self.extract_data(self)

    @property
    def variable_name(self):
        return self.variable_name

    @property
    def file_name(self):
        return self.file_name


# is it really useful?
def build_marker_default(probe):
    return "-"

# is it really useful?
def build_color_default(probe):
    return "red"


class DataView():
    def __init__(self):
        self.x_range = []
        self.y_range = []
        self.title = "empty title"
        self.x_label = "xLabel"
        self.y_label = "yLabel"
        self.probes = []

    def __repr__(self):
        repr_str = '''\txRange={}
        \tyRange={}
        \ttitle={}
        \txLabel={}
        \tyLabel={}
        '''.format(self.x_range, self.y_range, self.title, self.x_label, self.y_label)
        repr_str += "\tprobeList=\n"
        for p in self.probes:
            repr_str += "\t{}\n".format(p.name)
        return repr_str

# the following seems quite useless
# # TODO: wtf! use attribute or property
#     def set_xRange(self, xRange):
#         self.xRange = xRange

# # TODO: wtf! use attribute or property
#     def set_yRange(self, yRange):
#         self.yRange = yRange

# # TODO: wtf! use attribute or property
#     def set_xLabel(self, xLabel):
#         self.xLabel = xLabel

# # TODO: wtf! use attribute or property
#     def set_yLabel(self, yLabel):
#         self.yLabel = yLabel

# TODO: use join instead?

    def reset_probes(self):
        del self.probes
        self.probes = []

    def append_probes(self, new_probes):
        self.probes += new_probes

# TODO: wtf! use attribute or property
    def set_title(self, title):
        self.title = title

    def plot(self,
             view_settings={},
             build_label=build_label_default,
             build_marker=build_marker_default,
             build_color=build_color_default
             ):

        probe_view_settings = {k: {} for k in self.probes}
        logging.debug(probe_view_settings)

        master_tag_dict = {p.master_tag: p for p in self.probes}
        # TODO: can the following nested loop be recast into a single loop?
        for pattern, settings in view_settings.items():
            for p in dpath.search(master_tag_dict, "/" + pattern).values():
                probe_view_settings[p].update(settings)

        logging.debug(probe_view_settings)

        fig = plt.figure(figsize=(10, 10))
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.title)
        plt.grid(True)
        plt.gca().ticklabel_format(useOffset=False, style='sci')
        tickFormat = mtick.FormatStrFormatter("%.4e")
        plt.gca().get_yaxis().set_major_formatter(tickFormat)
        if len(self.x_range) > 0:
            plt.xlim(self.x_range[0], self.x_range[1])
        if len(self.y_range) > 0:
            plt.ylim(self.y_range[0], self.y_range[1])

        for probe in self.probes:
            logging.debug(probe.file_name
                          + " "
                          + probe.variable_name)
            x, var = probe.data
            # TODO : can we write this if-then in a more pythonist way?
            if build_label == None:
                probe_label = None
            else:
                probe_label = build_label(probe)
            if build_color == None:
                probe_color = None
            else:
                probe_color = build_color(probe)

            plot_args = probe_view_settings[probe]
#             plotArgs.update({'label' : probe_label})
            plot_args['label'] = probe_label
#            logging.debug(plotArgs)
            plt.plot(x, var, **plot_args)

        plt.legend(loc='best', frameon=False)
        plt.show(fig)


def plot_all_probes(db, variable_list, view_settings={}, x_range={}, y_range={}, x_label={}, y_label={}):
    if x_label == {}:
        x_label = {k: 'x' for k in variable_list}
    if y_label == {}:
        y_label = dict(zip(variable_list, variable_list))
    for variable in variable_list:
        probe_list = [f.get_default_probe()
                      for f in db if f.variable == variable]
        data_view = DataView()
        data_view.append_probes(probe_list)
        data_view.y_label = y_label.get(variable, [])
        data_view.x_label = x_label.get(variable, [])
        data_view.x_range = x_range.get(variable, [])
        data_view.y_range = y_range.get(variable, [])
        data_view.set_title(variable)
        data_view.plot(view_settings=view_settings)
