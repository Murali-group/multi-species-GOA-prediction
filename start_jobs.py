
# This script is meant to run scripts/jobs in 1 of 2 ways: 1) in an HPC environment or 2) screen jobs on local machine.
# It makes a config file for each algorithm and/or dataset being run.
# It has all the options of run_experiments.py, plus a some extra at the bottom. makes a config file for each

import yaml
import argparse
from collections import defaultdict
import os
import sys
import subprocess
import time
import itertools
import copy

# for some reason, the current directory isn't always in the path. this helps fix that
if "" not in sys.path:
    sys.path.insert(0,"")
from src.annotation_prediction.src import main as run_eval_algs
from src.annotation_prediction.src.algorithms import runner as runner
from src.utils import baobab_utils


def get_algs_to_run(alg_settings, **kwargs):
    # if there aren't any algs specified by the command line (i.e., kwargs),
    # then use whatever is in the config file
    if kwargs['algs'] is None:
        algs_to_run = run_eval_algs.get_algs_to_run(alg_settings)
        kwargs['algs'] = [a.lower() for a in algs_to_run]
        print("\nNo algs were specified. Using the algorithms in the yaml file:")
        print(str(kwargs['algs']) + '\n')
        if len(algs_to_run) == 0:
            print("ERROR: Must specify algs with --alg or by setting 'should_run' to [True] in the config file")
            sys.exit("Quitting")
    else:
        # make the alg names lower so capitalization won't make a difference
        kwargs['algs'] = [a.lower() for a in kwargs['algs']]
    return kwargs['algs']


def main(config_map, **kwargs):
    #input_settings = config_map['input_settings']
    #input_dir = input_settings['input_dir']
    algs = config_map['algs']
    #output_settings = config_map['output_settings']
    # combine the evaluation settings in the config file and the kwargs
    kwargs.update(config_map['eval_settings'])
    config_map = copy.deepcopy(config_map)
    #print("algs: %s" % (', '.join(a['name'] for a in algs)))
    kwargs['algs'] = get_algs_to_run(algs, **kwargs)

    yaml_file_pref = "" 
    if kwargs.get('split_taxon_prefix') is not None:
        # if split-taxons is specified, then put them in a subdir
        yaml_file_pref = "taxons/" 
        in_settings = config_map['input_settings']
        in_dir = in_settings['input_dir']
        stp = kwargs['split_taxon_prefix']
        for i in range(kwargs['split_taxon_ngroups'][0],kwargs['split_taxon_ngroups'][1]+1):
            taxon_file = "%s%s.txt" % (stp, i)
            print("Setting 'taxon_file' in the datasets to %s" % (taxon_file))
            if not os.path.isfile(taxon_file) and not os.path.isfile("%s/%s" % (in_dir, taxon_file)):
                sys.exit("ERROR: not found. Quitting")
            # the in_dir is added by the script, so remove that here
            taxon_file = taxon_file.replace("%s/"%in_dir,'') 
            for dataset in in_settings['datasets']:
                dataset['only_taxon_file'] = taxon_file
            curr_yaml_pref = "%s%s-" % (yaml_file_pref, i)
            run_jobs(algs, config_map, curr_yaml_pref, **kwargs)
    else:
        run_jobs(algs, config_map, yaml_file_pref, **kwargs)


def run_jobs(alg_settings, config_map, 
        yaml_file_pref='', postfix='', **kwargs):
    curr_config_map = copy.deepcopy(config_map)
    # setup the config file so that the specified algs have "True" for the should_run flag, and the others have false
    for alg, params in alg_settings.items():
        if alg.lower() in kwargs['algs']:
            #print('Running %s' % (alg))
            print(alg, params)
            params['should_run'] = [True]
        else:
            params['should_run'] = [False]
            continue
        # start one job per param combination
        if kwargs.get('job_per_param'):
            # get the parameter combinations
            combos = [dict(zip(params, val))
                for val in itertools.product(
                    *(params[param] for param in params))]
            for param_combo in combos:
                # only write the current alg, param settings in this yaml file
                curr_config_map['algs'] = {alg: {p: [val] for p, val in param_combo.items()}}
                # get the param str from the alg's runner 
                params_str = runner.get_runner_params_str(alg, {}, param_combo)
                curr_yaml_pref = yaml_file_pref+alg+params_str+postfix
                run_job_wrapper(curr_config_map, curr_yaml_pref, **kwargs)
    if not kwargs.get('job_per_param'):
        # run the specified algs together
        alg_name = '-'.join(kwargs['algs'])
        curr_yaml_pref = yaml_file_pref+alg_name+postfix
        run_job_wrapper(curr_config_map, curr_yaml_pref, **kwargs)


def run_job_wrapper(config_map, alg_name, **kwargs):
    # start a separate job for each dataset
    if kwargs.get('job_per_dataset'):
        for dataset in config_map['input_settings']['datasets']:
            curr_config_map = copy.deepcopy(config_map)
            # the exp_name is unique for each dataset.
            # just get the name and remove any folders(??)
            exp_name = dataset['exp_name'].split('/')[-1]
            net_version = dataset['net_version'].split('/')[-1]
            # add a folder for each dataset using the exp_name
            curr_alg_name = "%s-%s/%s" % (net_version, exp_name, alg_name)
            curr_config_map['input_settings']['datasets'] = [dataset]
            run_job(curr_config_map, curr_alg_name, **kwargs)
    else:
        run_job(config_map, alg_name, **kwargs)


def run_job(config_map, alg_name, **kwargs):
    # make a directory for this config, and then put a config file for each method inside
    yaml_base = kwargs['config'].replace('.yaml','')
    yaml_file = "%s/%s.yaml" % (yaml_base, alg_name)
    os.makedirs(os.path.dirname(yaml_file), exist_ok=True)
    cmd_file = os.path.abspath("%s/%s.sh" % (yaml_base, alg_name))
    log_file = os.path.abspath("%s/%s.log" % (yaml_base, alg_name))
    write_yaml_file(yaml_file, config_map)
    # now run it. Submit it to screen
    if kwargs.get('job_per_dataset'):
        name = alg_name
    else:
        name = "%s-%s" % (alg_name, config_map['input_settings']['datasets'][0]['exp_name'].split('/')[-1])
    # TODO add options for the python environment to use
    # default:
    #python = "/data/jeff-law/tools/anaconda3/bin/python"
    if kwargs.get('python'):
        python = kwargs['python']
    elif not kwargs.get('python') and kwargs.get('qsub') is True:
        # for timing/baobab:
        python = "source /data/jeff-law/projects/fungcat/2017-10-fungcat/virtual-envs/py3env-baobab/bin/activate \n" + \
                 "python"
    # for the csb machines:
    #python = "python"
    # pass the arguments specified when calling this script to this command
    str_args = get_script_args(**kwargs)
    command = "%s -u %s --config %s %s >> %s 2>&1" % (
        python, kwargs['script_to_run'], os.path.abspath(yaml_file), str_args, log_file)
    jobs = ["cd %s" % (os.getcwd()), command]
    if kwargs['qsub'] is True:
        submit = not kwargs['test_run']
        baobab_utils.writeQsubFile(
            jobs, cmd_file, name=name, submit=submit,  # log_file=log_file, # can't pass the log file since the PBS output file will overwrite 
            nodes=1, ppn=kwargs.get('cores',2), walltime=kwargs.get('walltime', '200:00:00'))
        # start a sleep job with the specified # cores and time to sleep
        if kwargs.get('sleep'):
            cores, sleep_time = kwargs['sleep']
            print("\tstarting sleep job with %s cores for %s time" % (cores, sleep_time))
            jobs = ["sleep %s" % sleep_time]
            sleep_file = "%s/sleep.sh" % (yaml_base)
            log_file = sleep_file.replace('.sh','.log')
            baobab_utils.writeQsubFile(
                jobs, sleep_file, name="sleep-%s"%sleep_time, submit=submit, log_file=log_file,
                nodes=1, ppn=int(cores), walltime=kwargs.get('walltime', '200:00:00'))
        if kwargs['test_run']:
            print(cmd_file)
            #sys.exit()
    else:
        # write the bash file
        write_bash_file(cmd_file, jobs)
        submit_to_screen(cmd_file, name, log_file, **kwargs)


# TODO make this more streamlined
def get_script_args(**kwargs):
    #str_args = "" 
    #for arg, val in kwargs.items():
    #    if arg in meta_args:
    #        continue
    #    str_args += " --%s %s" % (arg, val)
    #print(str_args)
    #sys.exit()
    if kwargs['pass_to_script']:
        args_to_pass = ' '.join(sys.argv[sys.argv.index('--pass-to-script')+1:])
    else:
        args_to_pass = ""
    return args_to_pass


def write_bash_file(cmd_file, jobs):
    print("\twriting to %s" % (cmd_file))
    with open(cmd_file, 'w') as out:
        out.write('echo "Job Started at: `date`"\n')
        # write each job, as well as an echo (print) statement of the job to be run to the qsub file
        out.write('\n'.join(['echo """%s"""\n%s' % (cmd, cmd) for cmd in jobs]) + '\n')
        out.write('echo "Job Ended at: `date`"\n')


def submit_to_screen(cmd_file, name, log_file, **kwargs):
    # looks like the job won't start if the name is too long, so truncate it here if needed
    name = name[:80] if len(name) > 80 else name
    # can't have '/' in screen name apparently 
    name = name.replace('/','-')
    print("\tsubmitting '%s' %s to screen" % (name, cmd_file))
    cmd = "screen -S %s -d -m /bin/sh -c \"bash %s >> %s 2>&1\"" % (name, cmd_file, log_file)
    print(cmd+'\n')
    if kwargs['test_run']:
        return
    else:
        subprocess.check_call(cmd, shell=True)


def write_yaml_file(yaml_file, config_map):
    print("\twriting to %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


# this is a list of the arguments that are unique to this script and should not be passed
# because the script we're running won't recognize them
meta_args = [
    "config",
    "script_to_run",
    "alg",
    "qsub",
    "test_run",
]


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description='Script for setting up various string experiments')

        parser.add_argument('--config', required=True,
            help='Configuration file')
        parser.add_argument('--taxon', '-T', dest="taxons", type=str, action='append',
                help="Specify the species taxonomy ID for which to evaluate. Multiple may be specified. Otherwise, all species will be used")
        #parser.add_argument('--string-

        parser.add_argument('--forced', action="store_true", default=False,
            help='Overwrite the ExpressionData.csv file if it already exists.')

    parser.add_argument('--skip-core', action='store_true', default=False,
            help="Skip running and evaluating the species in the core. " +
            "Used for the COMP and ELEC evaluations to run/eval only the species without EXPC annotations")

    parser.description = "This script is meant to run scripts/jobs in 1 of 2 ways: 1) in an HPC environment or 2) screen jobs on local machine. " + \
            "It makes a config file for each algorithm and/or dataset being run. " + \
            "It has all the options of run_experiments.py, plus a some extra at the bottom. makes a config file for each "

    group = parser.add_argument_group('start_jobs.py options')
    group.add_argument('--script-to-run', default="scripts/experiments.py",
            help="script to run when submitting to screen / qsub")
    group.add_argument('--alg', dest="algs", action="append", 
            help="Name of algorithm to run. May specify multiple. Default is whatever is set to true in the config file")
    group.add_argument('--job-per-param', action='store_true', default=False,
            help="Each parameter set combination per alg will get its own job")
    group.add_argument('--job-per-dataset', action='store_true', default=False,
            help="Each dataset will get its own job")
    group.add_argument('--split-taxon-prefix', type=str,
            help="Create a job per taxon file. This is the prefix of the taxon files. Should just be missing '#.txt'")
    group.add_argument('--split-taxon-ngroups', type=int, nargs=2,
            help="Start and end integers for the taxon groups to run (e.g., 1 10).")
    group.add_argument('--qsub', action='store_true', default=False,
            help="submit the jobs to a PBS queue with qsub. Currently only setup for GRISLI and SCINGE")
    group.add_argument('--test-run', action='store_true', default=False,
            help="Just print out the first command generated")
    group.add_argument('--cores', type=int, default=2,
            help="Number of cores to use per job submitted. Default: 2")
    group.add_argument('--sleep', type=str, nargs=2,
            help="<num-cores> <time-to-sleep> Amount of time to sleep (e.g., '12h') with the specified number of cores between submitted jobs. " + \
                    "Useful if jobs are RAM intensive at the start of the run, or timing methods on their own nodes")
    group.add_argument('--python', 
            help="Path to python bin executable to use. If qsub is specified, default is to use the baobab environment.")
    group.add_argument('--pass-to-script', action='store_true', default=False,
            help="All options specified after this option will be passed to the --script-to-run")
    # TODO add a machine option

    return parser


if __name__ == "__main__":
    # first load the run_eval_algs parser
    parser = run_eval_algs.setup_opts()
    parser = setup_parser(parser)
    opts = parser.parse_args()
    kwargs = vars(opts)
    kwargs['postfix'] = '' if kwargs['postfix'] is None else kwargs['postfix']
    config_file = opts.config
    with open(config_file, 'r') as conf:
        # for some reason this isn't recognized in other versions of PYyaml
        #config_map = yaml.load(conf, Loader=yaml.FullLoader)
        config_map = yaml.load(conf)
    main(config_map, **kwargs)
