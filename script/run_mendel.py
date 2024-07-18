import argparse
import math
import os
import subprocess
from dataclasses import asdict, dataclass
from pprint import pprint
import shutil
import multiprocessing
from functools import partial
import tqdm

from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pexpect
from jinja2 import Environment, FileSystemLoader, select_autoescape
from joblib import Memory
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 10:
    sys.exit('Python 3.10 or higher is required.')

my_dir = os.path.split(os.path.abspath(__file__))[0]

bin_dir = f'{my_dir}/../mendel/'
data_root_dir = f'{my_dir}/../data'
template_file = f'{bin_dir}/single_template.f'
cache_dir = './cache'

memory = Memory(cache_dir, verbose=0)

__all__ = [
    'write_combined_data', 'plot_relative_risk', 'plot_absolute_risk', 'run_mendel', 'bootstrap_mendel', 'memory',
    'bootstrap_run_mendel', 'get_lambda', 'base_frequencies'
]


def raise_helper(msg):
    raise Exception(msg)


base_frequencies = {
    'breast': [
        0.0220, 0.0140, 0.0050, 0.0000, 0.0010, 0.0020, 0.0030, 0.0040, 0.0040, 0.0030, 0.0020, 0.0010, 0.0020, 0.0180,
        0.0360, 0.0540, 0.0720, 0.1070, 0.2420, 0.3950, 0.5470, 0.6990, 0.9420, 1.7310, 2.6120, 3.4920, 4.3720, 5.4060,
        7.3600, 9.4670, 11.5750, 13.6820, 16.0160, 19.7130, 23.6360, 27.5580, 31.4810, 36.1100, 44.9770, 54.5500,
        64.1220, 73.6920, 83.0390, 91.0500, 98.8380, 106.6240, 114.4080, 121.7820, 126.7060, 131.2210, 135.7340,
        140.2440, 144.5710, 147.8110, 150.8680, 153.9230, 156.9750, 160.7430, 168.8230, 177.6160, 186.4040, 195.1860,
        204.0540, 213.4700, 222.9710, 232.4630, 241.9470, 250.8900, 256.6350, 261.8400, 267.0370, 272.2260, 276.4580,
        274.9920, 272.5720, 270.1490, 267.7210, 264.9980, 260.5260, 255.7610, 250.9930, 246.2230, 241.4770, 236.8880,
        232.3200, 227.7470, 223.1710, 218.1060, 210.1350, 201.6820
    ],
    'brain': [
        3.7450, 3.9430, 4.1700, 4.1870, 3.9950, 3.8020, 3.6100, 3.4240, 3.2740, 3.1300, 2.9860, 2.8420, 2.7030, 2.5980,
        2.4990, 2.4000, 2.3010, 2.2110, 2.1790, 2.1570, 2.1340, 2.1120, 2.1100, 2.2330, 2.3770, 2.5200, 2.6640, 2.8020,
        2.9140, 3.0200, 3.1270, 3.2330, 3.3330, 3.3910, 3.4430, 3.4940, 3.5460, 3.6050, 3.7070, 3.8170, 3.9260, 4.0360,
        4.1550, 4.3280, 4.5120, 4.6950, 4.8780, 5.0940, 5.5140, 5.9670, 6.4190, 6.8720, 7.3290, 7.8080, 8.2920, 8.7750,
        9.2580, 9.7450, 10.2540, 10.7670, 11.2790, 11.7920, 12.3220, 12.9590, 13.6130, 14.2670, 14.9210, 15.5690,
        16.1840, 16.7940, 17.4030, 18.0120, 18.5990, 19.0550, 19.4890, 19.9220, 20.3550, 20.7620, 21.0140, 21.2390,
        21.4640, 21.6890, 21.8390, 21.5450, 21.1770, 20.8100, 20.4420, 19.9830, 18.9740, 17.8740
    ],
    'colorectal': [
        0.0010, 0.0050, 0.0090, 0.0220, 0.0430, 0.0650, 0.0860, 0.1160, 0.1990, 0.2910, 0.3830, 0.4750, 0.5730, 0.7090,
        0.8510, 0.9930, 1.1350, 1.2730, 1.3820, 1.4870, 1.5920, 1.6960, 1.8240, 2.0860, 2.3720, 2.6570, 2.9430, 3.2620,
        3.7800, 4.3320, 4.8840, 5.4360, 6.0360, 6.9270, 7.8660, 8.8060, 9.7450, 10.7850, 12.4340, 14.1830, 15.9330,
        17.6820, 19.5540, 22.1620, 24.8920, 27.6220, 30.3510, 33.5010, 39.1800, 45.2780, 51.3740, 57.4690, 62.9320,
        64.6190, 65.6750, 66.7300, 67.7830, 69.2050, 72.8420, 76.8460, 80.8460, 84.8440, 89.0500, 94.5240, 100.2040,
        105.8790, 111.5490, 117.1350, 122.2510, 127.2830, 132.3070, 137.3240, 142.7000, 150.2750, 158.2030, 166.1160,
        174.0110, 181.8510, 189.4490, 196.9850, 204.4950, 211.9770, 219.2250, 225.2180, 230.9740, 236.6940, 242.3780,
        246.8160, 243.9830, 239.9320
    ],
    'leukemia': [
        4.5370, 6.1280, 7.9460, 8.3480, 7.3330, 6.3180, 5.3030, 4.3910, 4.1020, 3.9160, 3.7300, 3.5440, 3.3780, 3.3340,
        3.3100, 3.2860, 3.2620, 3.2320, 3.1650, 3.0920, 3.0190, 2.9460, 2.8860, 2.9050, 2.9370, 2.9690, 3.0010, 3.0420,
        3.1340, 3.2360, 3.3370, 3.4380, 3.5460, 3.6920, 3.8450, 3.9970, 4.1490, 4.3230, 4.6220, 4.9430, 5.2640, 5.5840,
        5.9250, 6.3860, 6.8660, 7.3470, 7.8280, 8.3480, 9.1060, 9.9040, 10.7020, 11.4990, 12.3370, 13.4160, 14.5360,
        15.6550, 16.7740, 17.9650, 19.5830, 21.2720, 22.9610, 24.6500, 26.4400, 28.8440, 31.3500, 33.8540, 36.3580,
        38.8850, 41.5580, 44.2550, 46.9490, 49.6420, 52.3860, 55.4460, 58.5560, 61.6630, 64.7660, 67.7850, 70.3220,
        72.7750, 75.2230, 77.6660, 80.0070, 81.7670, 83.4270, 85.0820, 86.7320, 88.0530, 87.4330, 86.4890
    ],
    'lfs': [
        11.8260, 12.7760, 13.9770, 13.9720, 12.7660, 11.5570, 10.3510, 9.2860, 9.0730, 9.0000, 8.9270, 8.8540, 8.7830,
        8.7280, 8.6750, 8.6220, 8.5700, 8.5150, 8.4460, 8.3760, 8.3050, 8.2350, 8.3130, 9.2950, 10.4270, 11.6570,
        12.7890, 14.0940, 16.4380, 18.9560, 21.4760, 23.9930, 26.7510, 30.9440, 35.3780, 39.8080, 44.2390, 49.4640,
        59.4410, 70.2110, 80.9800, 91.7460, 102.4250, 112.5790, 122.6450, 132.8100, 142.8720, 152.9390, 163.0570,
        173.1810, 183.3000, 193.4130, 203.9480, 217.0100, 230.4880, 244.0550, 257.5090, 271.9840, 292.6240, 314.2720,
        335.9980, 357.5990, 379.7300, 405.1800, 431.1490, 457.0740, 482.9540, 508.6620, 533.5790, 558.3060, 582.9520,
        607.5150, 630.8800, 647.4850, 662.8840, 678.1850, 693.3810, 706.6630, 709.0220, 709.5210, 709.9660, 710.3580,
        709.3740, 700.4610, 690.2430, 679.9390, 669.7530, 656.8650, 627.7350, 596.0340
    ]
}


@dataclass
class TemplateArgs:
    #
    age_cutoffs: tuple
    cancer_type: str

    # template_args can have
    travel: str = 'SEARCH'  # search or grid
    mxiter: int = 500
    #     number of age groups and model used (linear or piece wise)
    #
    # min parameter
    parmin: float = -5
    # max parameter
    parmax: float = 5
    # initial value
    parinit: float = 0.5

    maxage: int = 90
    #
    rr_model: str = 'linear'  # linear or piecewise

    npar: int = 0

    label: str = ''

    def __hash__(self):
        return hash((tuple(self.age_cutoffs), self.cancer_type, self.travel, self.mxiter, self.parmin, self.parmax,
                     self.parinit, self.rr_model))

    #
    def __post_init__(self) -> int:
        self.npar = len(self.age_cutoffs) - (0 if self.rr_model == 'linear' else 1)
        self.maxage = self.age_cutoffs[-1]
        assert self.maxage == 90, f'Last element of age_cutoff {self.age_cutoffs=} shold be 90'
        assert self.age_cutoffs[0] == 0, f'First element of age_cutoff {self.age_cutoffs=} should be 0'


def write_def_file(ped_name, ped_data):
    def_filename = os.path.join('mendel', ped_name, 'locusbr.dat')
    vars = list(set(ped_data['GT.gt.Variant_cDNA']))

    if not vars:
        raise ValueError(f'No variant for pedigree {ped_name}')

    # vars = vars[0]
    # def alleles(var):
    #     freq = 0.0001
    #     if var[-2] == '>':
    #         return freq, var[-3], var[-1]
    #     else:
    #         return freq, '_', 'D'

    with open(def_filename, 'w') as def_file:
        # for idx, var in enumerate(vars):
        #     if not isinstance(var, str):
        #         continue
        # freq, var1, var2 = alleles(var)
        freq = 0.0001
        def_file.write(f'''\
MAJOR   AUTOSOME 2 0
{1:<8}{1-freq:.6f}
{2:<8}{freq:.6f}
''')


def write_ped_file(dir_name, ped_data, exclude_proband, cancer_type):
    ped_data = ped_data.drop_duplicates(subset=['Pedigree name', 'UPN'])
    ped_filename = os.path.join(f'mendel/{dir_name}/pedigree.txt')

    locus_specific = False
    vars = set(ped_data['GT.gt.Variant_cDNA'])
    stats = defaultdict(int)

    if exclude_proband not in (None, False, 'all', 'first'):
        raise ValueError(f'Exclude proband can only be None, all or first')

    pedigrees = set(ped_data['Pedigree name'])
    print(f'Writing {ped_filename} for {len(pedigrees)} pedigree {dir_name} of size {ped_data.shape[0]}')
    with open(ped_filename, 'w') as ped_file:
        ped_file.write(f'''\
(i{len(str(ped_data.shape[0]))},a12)
(a7,4x,a7,2x,a7,2x,a2,4x,a2,4x,{6 + (1 if not locus_specific else len(vars))}(a4,2x),a4)
''')
        pd_idx = 0
        pd_name = None
        for idx, row in ped_data.iterrows():
            stats['# Peple'] += 1
            if row['Pedigree name'] != pd_name:
                pd_name = row['Pedigree name']
                pd_idx += 1
                pd_size = ped_data[ped_data['Pedigree name'] == pd_name].shape[0]
                ped_file.write(f'''{pd_size}{row['Pedigree name'][:8]:>12}\n''')
            # BASE
            age = 0 if np.isnan(row['Age']) else int(row['Age'])
            #
            # GENOTYPE
            #
            #gts = [''] * len(vars)
            if row['GT.Pos'] is True:
                gt = '1/2'
            elif row['GT.Neg'] is True:
                gt = '1/1'
            else:
                gt = ''

            # PHENOTYPE
            #

            if cancer_type == 'lfs':
                age1 = 0 if np.isnan(row['Breast_age_dx']) else int(row['Breast_age_dx'])
                age2 = 0 if np.isnan(row['Brain_age_dx']) else int(row['Brain_age_dx'])
                age3 = 0 if np.isnan(row['Sarcoma_age_dx']) else int(row['Sarcoma_age_dx'])
                age4 = 0 if np.isnan(row['Lung_age_dx']) else int(row['Lung_age_dx'])
                age5 = 0

                for cancer in ('Adrenal', 'Leukemia', 'Osteosarcoma', 'Phyllodes', 'Soft tissue sarcoma'):
                    if not np.isnan(row[cancer + '_age_dx']):
                        age5 = int(row[cancer + '_age_dx'])
                        break

                # age2 = 0 if np.isnan(row['Breast_age_dx']) else int(row['Breast_age_dx'])
                # age3 = 0 if np.isnan(row['Brain_age_dx']) else int(row['Brain_age_dx'])
                # age4 = 0 if np.isnan(row['Sarcoma_age_dx']) else int(row['Sarcoma_age_dx'])
                # age5 = 0
                # for cancer in ('Lung', 'Adrenal', 'Leukemia', 'Osteosarcoma', 'Phyllodes', 'Soft tissue sarcoma'):
                #     if not np.isnan(row[cancer + '_age_dx']):
                #         age5 = int(row[cancer + '_age_dx'])
                #         break
            else:
                age1 = 0 if np.isnan(row[f'{ cancer_type }_age_dx']) else int(row[f'{ cancer_type }_age_dx'])
                if age1 == 0 and row[cancer_type] is True:
                    age1 = 999
                # exclude proband for ascertainment correction
                if age1 != 0 and row['proband_flag_x'] == 'proband' and (exclude_proband == 'all' or \
                    (exclude_proband == 'first' and row[f'{ cancer_type }_age_dx'] <= row['Min.Dx.Age_Dx'])):
                    stats['proband_excluded'] += 1
                    age1 = 0

                # if cancer_type == 'Colorectal' and age1 < 18:
                #     age1 = 0
                age2 = 0
                age3 = 0
                age4 = 0
                age5 = 0

            agedeath = str(row['age at death'])
            if agedeath == 'nan':
                agedeath = '0'
            agelfu = int(age)
            # for cancer in ('Adrenal', 'Leukemia', 'Osteosarcoma', 'Phyllodes', 'Soft tissue sarcoma'):
            #      if not np.isnan(row[cancer + '_age_dx']) and row[cancer + '_age_dx'] > agelfu:
            #          agelfu = row[cancer + '_age_dx']

            fmt = '{upn:<7}    {mother:<7}  {father:<7}  {gender:<2}    {twin:<2}    {gt:<6}{age1:<4}  {age2:<4}  {age3:<4}  {age4:<4}  {age5:<4}  {agelfu:<4}  {agedeath:<4}\n'

            gender = row['Gender']
            if gender not in ('M', 'F'):
                # unknow, is she father or mother of someone?
                fathers = set(ped_data[ped_data['Pedigree name'] == row['Pedigree name']]['Father ID'])
                mothers = set(ped_data[ped_data['Pedigree name'] == row['Pedigree name']]['Mother ID'])
                if row['UPN'] in fathers:
                    gender = 'M'
                else:
                    gender = 'F'

            upn = row['UPN']
            mother = row['Mother ID'] if row['Mother ID'] else ''
            father = row['Father ID'] if row['Father ID'] else ''

            ped_file.write(
                fmt.format(
                    upn=str(upn),
                    mother=str(mother),
                    father=str(father),
                    gender=gender,
                    twin=' ',  # str(row['Twin status']) if row['Twin status'] > 0 else ' ',
                    gt=f'{gt.strip():<6}',
                    age1=age1,
                    age2=age2,
                    age3=age3,
                    age4=age4,
                    age5=age5,
                    agelfu=agelfu,
                    agedeath=agedeath,
                ))

    return stats


def write_combined_data(base_data, exclude_proband, name='combined', cancer_type='Colorectal'):
    # find all pedigrees
    dir = f'mendel/{name}'
    if not os.path.isdir(dir):
        os.mkdir(dir)

    stats = write_ped_file(name, base_data, exclude_proband=exclude_proband, cancer_type=cancer_type)
    write_def_file(name, base_data)
    return stats


def extract_results(result_file):
    sections = {
        0: [],
        1: [],
        2: [],
    }
    cur_section = 0

    with open(result_file) as summary:
        for line in summary:
            if line.strip().startswith('THE MAXIMUM'):
                cur_section = 1
            if line.strip().startswith('ASYMPTOTIC STANDARD ERRORS'):
                cur_section = 2
            if 'IERROR IGNORED' in line:
                continue
            sections[cur_section].append(line)
    #
    iter_line = sections[1][0]
    iter = iter_line.rstrip().rstrip(".").split()[-1]
    # print('\n'.join(x for x in sections[1] if x.strip()))
    print(f'return results from iteration {iter}')
    # look in section 0
    result = None
    for line in sections[0]:
        fields = line.strip().split()
        if len(fields) >= 4 and fields[0].isnumeric() and fields[1].isnumeric() and fields[0] == iter:
            result = [float(x.replace('D+', 'E').replace('D-', 'E-').replace('D', 'E')) for x in line.strip().split()
                     ][3:]
            # print(f'START {line=}, {result=}')
            continue
        if result:
            if not line.strip():
                break
            result.extend(
                [float(x.replace('D+', 'E').replace('D-', 'E-').replace('D', 'E')) for x in line.strip().split()])
            # print(f'ADD {line=}')
    #
    sd = None
    for line in sections[2]:
        fields = line.strip().split()
        # 0.3199D-01
        if sd:
            if not line.strip():
                break
            sd.extend([
                np.nan if x == 'NaN' else float(x.replace('D+', 'E').replace('D-', 'E-').replace('D', 'E'))
                for x in line.strip().split()
            ])
            #print(f'ADD SD {line=}')
            continue
        if len(fields) > 0 and (fields[0][0].isnumeric() or fields[0] == "NaN"):
            sd = [
                np.nan if x == "NaN" else float(x.replace('D+', 'E').replace('D-', 'E-').replace('D', 'E'))
                for x in line.strip().split()
            ]
            #print(f'SD {line=}')
            continue

    # result = prepend + result + postpend

    if not sd or all(np.isnan(sd)):
        sd = []

    if result:
        return {'rr': result, 'sd': sd}

    print(f'No result is found from {result_file}')
    return {'rr': [], 'sd': []}


def generate_source_code(template_args):
    if not os.path.isfile(template_file):
        raise RuntimeError(f'No template {template_file}')

    env = Environment(loader=FileSystemLoader(os.path.split(template_file)[0]), autoescape=select_autoescape())
    env.globals['raise'] = raise_helper

    template = env.get_template(os.path.basename(template_file))
    source_code = template.render(**asdict(template_args))
    prog_name = f'single{hash(template_args)}'

    if not os.path.isfile(f'{bin_dir}/{prog_name}'):
        with open(f'{bin_dir}/{prog_name}.f', 'w') as src:
            src.write(source_code)

        compile_cmd = f'gfortran -O3 mendela.f {prog_name}.f -o {prog_name}'

        print(f'Compiling program {prog_name}')
        subprocess.run(
            compile_cmd, cwd=bin_dir, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return f'{bin_dir}/{prog_name}'


@memory.cache
def run_mendel(data_dir, keep_results=True, **kwargs):
    #
    bootstrap_pedigrees = kwargs.pop('bootstrap_pedigrees', False)
    template_args = TemplateArgs(**kwargs)

    if not os.path.isdir(data_dir):
        data_dir = os.path.join(data_root_dir, data_dir)

    if not os.path.isdir(data_dir):
        raise ValueError(f'Non existing data dir {data_dir}')

    ped_name = os.path.basename(data_dir)
    cmd = generate_source_code(template_args)

    tdir = f'{data_dir}{hash(template_args)}'
    if bootstrap_pedigrees and os.path.isdir(tdir):
        # another process has processed this particular argument
        raise RuntimeError('Ignore bootstrap configuration.')

    print(f'Processing {ped_name} at {tdir}')
    if not os.path.isdir(tdir):
        os.makedirs(tdir)
        os.symlink(os.path.join(data_dir, 'locusbr.dat'), os.path.join(tdir, 'locusbr.dat'))
        # get number of pedigrees
        if bootstrap_pedigrees:
            with open(os.path.join(data_dir, 'pedigree.txt')) as ifile, \
                open(os.path.join(tdir, 'pedigree.txt'), 'w') as ofile:
                # first two lines
                ofile.write(ifile.readline())
                ofile.write(ifile.readline())
                pedigrees = []
                for line in ifile:
                    if len(line.split()) == 2:
                        # try to avoid duplicated name
                        nfam = line.split()[0]
                        pedigrees.append([f'{int(nfam): <5}PED{len(pedigrees)+ 1}\n'])
                        continue
                    pedigrees[-1].append(line)
                # write pedigrees.
                choices = np.random.choice(len(pedigrees), size=len(pedigrees), replace=True)
                choices.sort()
                for choice in choices:
                    ofile.write(''.join(pedigrees[choice]))
        else:
            os.symlink(os.path.join(data_dir, 'pedigree.txt'), os.path.join(tdir, 'pedigree.txt'))

    result_file = f'{tdir}/single_out.dat'
    proc = pexpect.spawn([cmd], cwd=tdir, timeout=500)

    # DO YOU WISH TO CHANGE TO BATCH MODE? [YES/NO]
    proc.expect(':')
    proc.sendline('no')
    # CHOOSE AN ITEM [1,...,21]:
    proc.expect(':')
    proc.sendline('21')
    # ANOTHER PROBLEM [YES/NO]:
    proc.expect(':')
    proc.sendline('no')
    # proc.interact()
    proc.expect(pexpect.EOF)

    # os.remove(cmd[:-2])
    # os.remove(cmd)
    res = extract_results(result_file) | asdict(template_args)
    # print(f'random samples kept in {tdir}')
    if not keep_results:
        shutil.rmtree(tdir)
    return res


def bootstrap_run_mendel(idx, data_dir, **kwargs):
    #
    if not os.path.isdir(data_dir):
        data_dir = os.path.join(data_root_dir, data_dir)

    if not os.path.isdir(data_dir):
        raise ValueError(f'Non existing data dir {data_dir}')

    # n_age_groups = len(age_cutoffs) - 1
    max_age = 90
    np.random.seed(idx)

    while True:
        raw_ags = np.random.randint(10, 30, 100)
        age_cutoffs = raw_ags.cumsum()
        age_cutoffs = [0] + list(age_cutoffs[:sum(age_cutoffs < max_age)])
        if max_age - age_cutoffs[-1] > 5:
            age_cutoffs.append(max_age)
        else:
            age_cutoffs[-1] = max_age

        if len(age_cutoffs) <= 7:
            break

    kwargs['age_cutoffs'] = age_cutoffs
    kwargs['bootstrap_pedigrees'] = True
    kwargs['label'] = ''
    try:
        return run_mendel(data_dir, **kwargs) | {'bootstrap_age_cutoffs': age_cutoffs, 'bootstrap_index': idx}
    except (FileNotFoundError, RuntimeError):
        # sometimes race condition causes two instances with the same age groups compete with
        # each other and one of them remove the files after done.
        return None


@memory.cache
def bootstrap_mendel(data_dir, bootstrap, njobs=5, **kwargs):
    results = []

    with multiprocessing.Pool(njobs) as pool:
        starting = 0
        complted_age_cutoffs = set()
        while True:
            n_missing = bootstrap - len(results)
            starting += n_missing
            if n_missing == 0:
                break
            # for res in list(
            #         tqdm.tqdm(
            # total=n_missing)
            for res in pool.map(
                    partial(bootstrap_run_mendel, data_dir=data_dir, **kwargs), range(starting, n_missing + starting)):
                if not res:
                    # failed, possibly due to overlap
                    continue
                if not res['sd']:
                    print(f'Ignore result with age group {res["age_cutoffs"]} due to missing SD')
                    continue
                if any(x == 0 for x in res['sd']) or any(x in (res['parmin'], res['parmax']) for x in res['rr']):
                    print(f'Ignore result hitting bounday: {res["age_cutoffs"]} {res["rr"]}')
                    continue
                if tuple(res['age_cutoffs']) in complted_age_cutoffs:
                    # print('duplicated age group')
                    continue
                else:
                    results.append(res)
                    complted_age_cutoffs.add(tuple(res['age_cutoffs']))
                    print(f'{len(results)}/{bootstrap}: {res["age_cutoffs"]}. {res["rr"]}')

    # now, let us summarize
    all_grid = None
    for result in results:
        xx = range(result['maxage'])
        indexes = [[idx for idx, x in enumerate(result['age_cutoffs'][1:]) if age < x][0] for age in xx]
        if result['rr_model'] == 'linear':
            age_cutoff = result['age_cutoffs']
            y = [
                result['rr'][idx] + (x - age_cutoff[idx]) * (result['rr'][idx + 1] - result['rr'][idx]) /
                (age_cutoff[idx + 1] - age_cutoff[idx]) for x, idx in zip(xx, indexes)
            ]
        else:
            y = [result['rr'][idx] for idx in indexes]

        if all_grid is None:
            all_grid = np.array([y])
        else:
            all_grid = np.vstack([all_grid, y])

    # calculate mean and confidence interval
    return {
        'age_cutoffs': range(results[0]['maxage'] + 1),
        'cancer_type': results[0]['cancer_type'],
        'rr_model': 'bootstrap',
        'rr': all_grid.mean(axis=0),
        'rr_low': list(pd.DataFrame(all_grid).quantile(axis=0, q=0.025)),
        'rr_high': list(pd.DataFrame(all_grid).quantile(axis=0, q=0.975)),
        'bootstrap': bootstrap,
        'max_age': results[0]['maxage'],
        'raw_results': results,
        'label': kwargs.get('label', '')
    }


def plot_relative_risk(results,
                       title='',
                       ymin=None,
                       ymax=None,
                       ylog=False,
                       xlabel=None,
                       ylabel=None,
                       legends=None,
                       plotter=None):
    if not plotter:
        plotter = plt

    plt.style.use('ggplot')
    maxy = 0
    miny = 100
    for idx, result in enumerate(results):
        xx = range(result['age_cutoffs'][-1])
        if result['rr_model'] != 'bootstrap' and not result['sd']:
            result['sd'] = [0] * len(result['rr'])
        indexes = [[idx for idx, x in enumerate(result['age_cutoffs'][1:]) if age < x][0] for age in xx]
        if result['rr_model'] == 'linear':
            age_cutoff = result['age_cutoffs']
            y = [
                math.exp(result['rr'][idx] + (x - age_cutoff[idx]) * (result['rr'][idx + 1] - result['rr'][idx]) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]
            y_low = [
                math.exp(result['rr'][idx] - (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96) +
                         (x - age_cutoff[idx]) *
                         (result['rr'][idx + 1] -
                          (0 if np.isnan(result['sd'][idx + 1]) else result['sd'][idx + 1] * 1.96) - result['rr'][idx] +
                          (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96)) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]
            y_high = [
                math.exp(result['rr'][idx] + (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96) +
                         (x - age_cutoff[idx]) *
                         (result['rr'][idx + 1] +
                          (0 if np.isnan(result['sd'][idx + 1]) else result['sd'][idx + 1] * 1.96) - result['rr'][idx] -
                          (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96)) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]

        elif result['rr_model'] == 'piecewise':
            y = [math.exp(result['rr'][idx]) for idx in indexes]
            y_low = [
                math.exp(result['rr'][idx] - (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96))
                for idx in indexes
            ]
            y_high = [
                math.exp(result['rr'][idx] + (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96))
                for idx in indexes
            ]
        else:
            y = [math.exp(x) for x in result['rr']]
            y_low = [math.exp(x) for x in result['rr_low']]
            y_high = [math.exp(x) for x in result['rr_high']]

        if legends and len(legends) > idx:
            label = legends[idx]
        else:
            label = result.get('label', '') or ' - '.join(str(x) for x in result['age_cutoffs'])

        plotter.plot(xx, y, label=label)
        plotter.fill_between(xx, y_low, y_high, color='b', alpha=.1)

        if len(result['age_cutoffs']) == 2:
            plotter.text(60, y[0] + 2, f'rr={y[0]:.1f} (CI: {y_low[0]:.1f} - {y_high[0]:.1f})')
        maxy = max(maxy, max(y_high))
        miny = min(miny, min(y_low))

    # plotter.title(title)
    plotter.hlines(1, xmin=0, xmax=90, colors='green', linestyles='dotted')
    # plotter.text(72, 1.5, 'Relative risk = 1')
    plotter.legend()
    if plotter == plt:
        plotter.title(title)
        plotter.xlabel(xlabel or 'Relative Risk')
        plotter.ylabel(ylabel or 'Age at Diagnosis')
        plotter.ylim(ymin if ymin else miny * 0.9, ymax if ymax else maxy * 1.1)
        if ylog:
            plotter.yscale("log")
        else:
            plotter.yticks(list(plotter.yticks()[0]) + [1])
    else:
        plotter.set_title(title)
        plotter.set_xlabel(xlabel or 'Relative Risk')
        plotter.set_ylabel(ylabel or 'Age at Diagnosis')
        plotter.set_ylim(ymin if ymin else miny * 0.9, ymax if ymax else maxy * 1.1)
        if ylog:
            plotter.set_yscale("log")
        else:
            plotter.set_yticks(list(plotter.get_yticks()) + [1])
    return plotter


def get_lambda(rrbr, popbr):
    p1 = (1 - 0.0004)**2
    p2 = 1 - p1

    ffncbr = [1] + [0] * 90
    ffbr = [1] + [0] * 90

    cumncbr = 0
    cumbrrisk = 0
    lambda0 = [0] * 90
    lambda1 = [0] * 90

    for iage in range(90):
        lambda0[iage]= (popbr[iage]/100000.0)*(p1*ffncbr[iage] + p2*ffbr[iage])/ \
                        (p1*ffncbr[iage] + p2*rrbr[iage]*ffbr[iage])
        lambda1[iage] = lambda0[iage] * rrbr[iage]
        cumncbr = cumncbr + lambda0[iage]
        ffncbr[iage + 1] = math.exp(-cumncbr)
        cumbrrisk = cumbrrisk + lambda1[iage]
        ffbr[iage + 1] = math.exp(-cumbrrisk)
    return [100000 * x for x in lambda0], [x * 100000 for x in lambda1]


def plot_absolute_risk(results, title='', ymin=None, ymax=None, legends=None, xlabel=None, ylabel=None, plotter=None):
    popbr = base_frequencies[results[0]['cancer_type']]

    if not plotter:
        plotter = plt

    # plt.plot(popbr, label='Population incidence rate')
    maxy = 0
    miny = 100
    plt.style.use('ggplot')
    for result in results:
        xx = range(len(popbr))
        indexes = [[idx for idx, x in enumerate(result['age_cutoffs'][1:]) if age < x][0] for age in xx]
        if result['rr_model'] == 'linear':
            age_cutoff = result['age_cutoffs']
            y = [
                math.exp(result['rr'][idx] + (x - age_cutoff[idx]) * (result['rr'][idx + 1] - result['rr'][idx]) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]
            y_low = [
                math.exp(result['rr'][idx] - (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96) +
                         (x - age_cutoff[idx]) *
                         (result['rr'][idx + 1] -
                          (0 if np.isnan(result['sd'][idx + 1]) else result['sd'][idx + 1] * 1.96) - result['rr'][idx] +
                          (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96)) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]
            y_high = [
                math.exp(result['rr'][idx] + (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96) +
                         (x - age_cutoff[idx]) *
                         (result['rr'][idx + 1] +
                          (0 if np.isnan(result['sd'][idx + 1]) else result['sd'][idx + 1] * 1.96) - result['rr'][idx] -
                          (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96)) /
                         (age_cutoff[idx + 1] - age_cutoff[idx])) for x, idx in zip(xx, indexes)
            ]
        elif result['rr_model'] == 'piecewise':
            y = [math.exp(result['rr'][idx]) for idx in indexes]
            y_low = [
                math.exp(result['rr'][idx] - (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96))
                for idx in indexes
            ]
            y_high = [
                math.exp(result['rr'][idx] + (0 if np.isnan(result['sd'][idx]) else result['sd'][idx] * 1.96))
                for idx in indexes
            ]
        else:
            y = [math.exp(x) for x in result['rr']]
            y_low = [math.exp(x) for x in result['rr_low']]
            y_high = [math.exp(x) for x in result['rr_high']]
        #

        lambda0, lambda1 = get_lambda(y, popbr)
        lambda0_high, lambda1_high = get_lambda(y_high, popbr)
        lambda0_low, lambda1_low = get_lambda(y_low, popbr)

        maxy = max(maxy, max(lambda1_high))
        miny = min(miny, min(lambda1_low))

        if legends is None:
            legends = {}
        p0 = plotter.plot(
            xx,
            lambda0,
            label=legends.get('non-carriers', 'Non-carriers: ' + result.get('label', 'TP53 mutation carriers')))
        p1 = plotter.plot(
            xx,
            lambda1,
            label=legends.get('carriers', 'TP53 Carriers: ' + result.get('label', 'TP53 mutation carriers')))
        plotter.fill_between(xx, lambda0_low, lambda0_high, color=p0[0].get_color(), alpha=.1)
        plotter.fill_between(xx, lambda1_low, lambda1_high, color=p1[0].get_color(), alpha=.1)
        # maxy = max(maxy, max(lambda1))

    plotter.legend()
    if plotter == plt:
        plotter.title(title)
        plotter.ylabel(ylabel or 'Rate per 100,000')
        plotter.xlabel(xlabel or 'Age at Diagnosis')
        plotter.ylim(ymin if ymin else miny * 0.9, ymax if ymax else maxy * 1.1)
    else:
        plotter.set_title(title)
        plotter.set_ylabel(ylabel or 'Rate per 100,000')
        plotter.set_xlabel(xlabel or 'Age at Diagnosis')
        plotter.set_ylim(ymin if ymin else miny * 0.9, ymax if ymax else maxy * 1.1)

    return plotter


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--template', default=template_file, help='full path to the template code')
    parser.add_argument(
        '--rr-model', default='linear', choices=['linear', 'piecewise'], help='Model for the relative risk.')
    parser.add_argument('--travel', default='SEARCH', choices=['SEARCH', 'GRID'], help='''Search method''')
    parser.add_argument('--mxiter', default=500, type=int, help='''Max iterations''')
    parser.add_argument('--parinit', default=0.5, help='''Initial value for parameter''')
    parser.add_argument('--parmin', default=0, help='''Minimal value for parameter''')
    parser.add_argument('--parmax', default=5, help='''Max parameter''')
    parser.add_argument(
        '--cancer-type',
        required=True,
        choices=['leukemia', 'colorectal', 'breast', 'brain', 'lfs'],
        help='''Cancer type''')
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        '--age-cutoffs',
        type=int,
        nargs='+',
        help='''Cutoff values, the first value should be zero. The last value
            should be is the maximum age, and patients with age greater
            than the range will be censored.''')
    action.add_argument(
        '--bootstrap',
        type=int,
        help='''If set to True, will bootstrap the cutoff ages. The first and last
            cutoff values will be fixed, the middle values will be randomized for specified
            times.''')
    parser.add_argument('--data-dir', required=True, help='Data directory')
    parser.add_argument('--label', help='Comment for the analysis.')
    args = parser.parse_args()

    if args.age_cutoffs is None:
        res = bootstrap_mendel(
            args.data_dir,
            travel=args.travel,
            mxiter=args.mxiter,
            parmin=args.parmin,
            parmax=args.parmax,
            parinit=args.parinit,
            age_cutoffs=args.age_cutoffs,
            cancer_type=args.cancer_type,
            rr_model=args.rr_model,
            label=args.label or '',
            bootstrap=args.bootstrap)
    else:
        res = run_mendel(
            args.data_dir,
            travel=args.travel,
            mxiter=args.mxiter,
            parmin=args.parmin,
            parmax=args.parmax,
            parinit=args.parinit,
            age_cutoffs=args.age_cutoffs,
            cancer_type=args.cancer_type,
            rr_model=args.rr_model,
            label=args.label or '',
        )
    pprint(res)
