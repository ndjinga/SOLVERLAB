import resultsBase as rb
import re
import numpy as np
import VTU_routines
import MED_routines
import CoreFlows_variables_dictionary

solver_dict_cf = {
    'cf-ordre1': 'CF-ordre1',
    'vfroe_': 'CF-ordre1',
    'vfroecorp': 'CF-ordre2',
    'cf-ordre2': 'CF-ordre2',
    'vffc_wb': 'VFFCeq',
    'pstag_wb': 'pSTAGeq',
    'ps_wb': 'pSTAGeq',
    'staggered_wb': 'pSTAGeq'
}

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


def get_variable_name_from_string(str, varDict, default_value="noVar"):
    res = default_value
    for var_name in varDict.keys():
        if re.search(var_name, str,re.IGNORECASE):
            res = varDict[var_name]
            break
    return res


# deprecated: use the last_step property instead
def get_last_step(resDb):
    res = max([elt.step for elt in resDb])
    return res

# deprecated: use last_step_elts property instead
def extract_last_step_elt(resDb):
    last_step = get_last_step(resDb)
    last_step_elts = [elt for elt in resDb if elt.step == last_step]
    if len([elt for elt in last_step_elts if elt.is_stationnary == True]) > 0:
        logging.warning(
            "you are extracting last step variables "
            "from a result list that contains stationnary results.\n"
        )
    return last_step_elts


def is_stationnary_cf(name):
    pattern_for_stationnary = ".*_Stat_.*"
    match = re.search(pattern_for_stationnary, name)
    return (match != None)


def create_result_list_1D(dir_list, tags=[], fname_glob=None, code='cf', solver=None, search_step=True):
    if len(dir_list) == 0:
        logging.warning(
            "You passed an empty list of dir to CreateResultList_CoreFlows")
        return []
    else:
        result_list = []
        for a_dir in dir_list:
            if fname_glob == None:
                fileList = rb.create_file_list(a_dir, ext='csv')
            else:
                fileList = rb.create_file_list(a_dir, fname_glob=fname_glob)
            for f in fileList:
                variable = rb.get_variable_name_from_string(
                    f, var_dict_cf, default_value='unknown_variable')
                if search_step :
                    step = rb.get_step_from_root_name(f)
                else:
                    step = None
                is_stationnary = is_stationnary_cf(f)
                if solver == None:
                    result_solver = rb.get_solver_name(f, solver_dict_cf)
                else:
                    result_solver = solver
                field_result = rb.Result(
                    file_name=f,
                    variable=variable,
                    dimension=1,
                    step=step,
                    code=code,
                    is_stationnary=is_stationnary,
                    solver=result_solver,
                    tags=tags
                )
                result_list.append(field_result)
    return result_list


def create_result_list_2D(dir_list, tags=[], solver=None):
    if len(dir_list) == 0:
        logging.warning(
            "You passed an empty list of dir to CreateResultList_CoreFlows")
        return []
    else:
        result_list = []
        for a_dir in dir_list:
            files = rb.create_file_list(a_dir, ext='vtu')
            for f in files:
                variable = rb.get_variable_name_from_string(
                    f, var_dict_cf_2D, default_value='unknown_variable')
                step = rb.get_step_from_root_name(f)
                is_stationnary = is_stationnary_cf(f)
                if solver == None:
                    result_solver = rb.get_solver_name(f, solver_dict_cf)
                else:
                    result_solver = solver
                field_result = rb.Result(
                    file_name=f,
                    variable=variable,
                    dimension=2,
                    step=step,
                    code='cf',
                    is_stationnary=is_stationnary,
                    solver=result_solver,
                    tags=tags
                )
                result_list.append(field_result)
    return result_list

# for the sake of compatibility with ancient version
CreateResultList_2D = create_result_list_2D
