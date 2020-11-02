import copy
import itertools as it
import json
import os
from optparse import OptionParser, Values
from socket import gethostname
from typing import Dict, List


def get_option(options: Values, key: str):
    """Getting paths and the name of the configuration files, and the number of cores.

    :param key: argument name corresponding to path or file name / number of cores
    :param options: contains the value of the argument
    """
    if not hasattr(options, key):
        raise ValueError(f"There is no {key} option")
    return getattr(options, key)


def get_git_tag(path: str) -> str:
    prev_wd = os.getcwd()
    os.chdir(path)

    command = "git rev-parse --short HEAD"
    git_hash = os.popen(command).read().strip()

    os.chdir(prev_wd)
    return git_hash


def assign_level(nested_dicts, comp_key, val):
    d = nested_dicts
    for k in comp_key[:-1]:
        d = d[k]
    d[comp_key[-1]] = val


def apply_patch(conf: Dict, patch: Dict[str, str]):
    for key, val in patch.items():
        tuple_key = tuple(key.split("->"))
        assign_level(conf, tuple_key, val)


def create_configs(
    task_filename: str,
    configs_folder: str,
    parent_config_filename: str,
    git_tag: str,
    counter_file: str,
    results_folder: str,
    machine_id: str,
):
    """Creates parent_config files based on the parent parent_config file and writes to them all possible
    combinations of parameter values from the job file. .

    :return: a list with the names of the created files
    """
    list_configs = []
    with open(task_filename) as f:
        task = json.load(f)
    with open(parent_config_filename) as f:
        parent_config = json.load(f)

    count = get_starting_number_from_file(counter_file) - 1
    cross_scan_gen = (cross_scan_dict for cross_scan_dict in task["Task"]["CrossScan"] if cross_scan_dict)
    for cross_scan_dict in cross_scan_gen:
        patches = (dict(zip(cross_scan_dict.keys(), comb)) for comb in it.product(*cross_scan_dict.values()))
        for patch in patches:
            count += 1
            config_filename = f"{count}.json"
            save_folder = os.path.join(results_folder, f"{count}_{machine_id}/")
            if not os.path.exists(save_folder):
                os.mkdir(save_folder)

            config = copy.deepcopy(parent_config)
            apply_patch(config["Configuration"], patch)
            apply_patch(
                config,
                {
                    "Configuration->gitHashShort": git_tag,
                    "Configuration->Name": str(count),
                    "Configuration->LoggerParameters->FilePath": save_folder,
                    "Configuration->ParentConfigName": parent_config_filename,
                    "Configuration->TaskFile": task_filename,
                    "Configuration->ComputerID": machine_id,
                },
            )
            with open(os.path.join(configs_folder, config_filename), "w") as f:
                json.dump(config, f, sort_keys=True, indent=2)
            list_configs.append(config_filename)

    set_counter_file(counter_file, count + 1)
    return list_configs


def create_params_file(filename: str, config_files_list: List[str]) -> None:
    with open(filename, "w") as f:
        f.writelines(map(lambda x: f"{x}\n", config_files_list))


def get_starting_number_from_file(filename: str) -> int:
    with open(filename) as f:
        lines = f.readlines()
    return int(lines[0].strip())


def set_counter_file(filename: str, count: int) -> None:
    with open(filename, "w") as f:
        f.write(f"{count}\n")


def get_machine_id() -> str:
    return gethostname()


def main():
    parser = OptionParser()

    parser.add_option("-t", "--task_file_name")
    parser.add_option("-p", "--parent_config_filename")
    parser.add_option("-r", "--results_folder")
    parser.add_option("-g", "--git_path")
    parser.add_option("-c", "--counter_filename")
    parser.add_option("-f", "--configs_folder")

    (options, args) = parser.parse_args()

    results_folder = get_option(options, "results_folder")
    parent_config_filename = get_option(options, "parent_config_filename")
    task_filename = get_option(options, "task_file_name")
    git_path = get_option(options, "git_path")
    counter_filename = get_option(options, "counter_filename")
    configs_folder = get_option(options, "configs_folder")

    machine_id = get_machine_id()

    config_files = create_configs(
        task_filename=task_filename,
        configs_folder=configs_folder,
        parent_config_filename=parent_config_filename,
        git_tag=get_git_tag(git_path),
        counter_file=counter_filename,
        results_folder=results_folder,
        machine_id=machine_id,
    )

    create_params_file(os.path.join(configs_folder, "sims.dat"), config_files)
    print("\n".join(config_files))


if __name__ == "__main__":
    main()
