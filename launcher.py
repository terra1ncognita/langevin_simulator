import copy
import itertools as it
import json
import os
from optparse import OptionParser, Values
from socket import gethostname
from typing import Any, Dict, List

import numpy as np
from scipy.optimize import fsolve


class MTSpring:
    def __init__(self, kw, wb, a, b, c, sb, ks, os):
        self.kw = kw
        self.wb = wb
        self.a = a
        self.b = b
        self.c = c
        self.sb = sb
        self.ks = ks
        self.os = os

    def __call__(self, x: float) -> float:
        if x < 0:
            return 0
        elif 0 <= x < self.wb:
            return self.kw * x
        elif self.wb <= x < self.sb:
            return self.a * x ** 2 + self.b * x + self.c
        else:
            return self.ks * x + self.os

    @classmethod
    def from_config(cls, config: Dict[str, Any], side: str) -> "MTString":
        mp = config['Configuration']["ModelParameters"]
        if side == "R":
            return cls(
                kw=mp["MTstiffWeakSlopeL"],
                wb=mp["MTstiffWeakBoundaryL"],
                a=mp["MTstiffParabolicAL"],
                b=mp["MTstiffParabolicBL"],
                c=mp["MTstiffParabolicCL"],
                sb=mp["MTstiffStrongBoundaryL"],
                ks=mp["MTstiffStrongSlopeL"],
                os=mp["MTstiffStrongIntersectL"],
            )
        elif side == "L":
            return cls(
                kw=mp["MTstiffWeakSlopeR"],
                wb=mp["MTstiffWeakBoundaryR"],
                a=mp["MTstiffParabolicAR"],
                b=mp["MTstiffParabolicBR"],
                c=mp["MTstiffParabolicCR"],
                sb=mp["MTstiffStrongBoundaryR"],
                ks=mp["MTstiffStrongSlopeR"],
                os=mp["MTstiffStrongIntersectR"],
            )
        else:
            raise ValueError(f"side must be either R or L, got {side}")


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


def assign_initial_conditions(config: Dict[str, Any]) -> None:
    """Compute initial values for all coordinated and assign to config inplace."""
    f, f_pre, k_r, k_l, mt_lenght = [
        config['Configuration']["ModelParameters"][x]
        for x in ["movementTotalForce", "prestretchTotalForce", "trapstiffR", "trapstiffL", "MTlength"]
    ]

    mt_spring_r = MTSpring.from_config(config, side="R")
    mt_spring_l = MTSpring.from_config(config, side="L")

    def prestretch_calc(y):
        d1, d2 = y
        return [mt_spring_l(2 * d1) - f_pre, mt_spring_r(2 * d2) - f_pre]

    # d1 > 0, d2 > 0
    d1, d2 = fsolve(prestretch_calc, x0=[0, 0])

    x_bead_l = -d1 - mt_lenght / 2.0
    x_bead_r = d2 + mt_lenght / 2.0
    x_trap_l = -f_pre / k_l + x_bead_l
    x_trap_r = f_pre / k_r + x_bead_r

    apply_patch(
        config,
        {
            "Configuration->InitialConditions->direction": "1.0",
            "Configuration->InitialConditions->xPed": "0.0",
            "Configuration->InitialConditions->xMol": "0.0",
            "Configuration->InitialConditions->xMT": "0.0",
            "Configuration->InitialConditions->xBeadl": str(x_bead_l),
            "Configuration->InitialConditions->xBeadr": str(x_bead_r),
            "Configuration->InitialConditions->xTrapl": str(x_trap_l),
            "Configuration->InitialConditions->xTrapr": str(x_trap_r),
            "Configuration->InitialConditions->phi": str(np.pi / 3),
        },
    )


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
            assign_initial_conditions(config)

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
