#!/usr/bin/python3

import os, sys, re, json
from collections import defaultdict as dd, OrderedDict as od

POP = ['EUR', 'AFR', 'FIN', 'ARAB', 'SAS', 'EAS', 'OTH', 'CSA', 'AMR', 'NON-EUR', 'HIS', 'ALL']
ANA = ['A1', 'A2', 'B1', 'B2', 'B3', 'C1', 'C2', 'D1']
AGE = ['ALL', 'GT_60', 'LE_60']

def parse_file(uri):
    filename = uri.split('/').pop()
    filename = re.sub('.txt|.gz|.bgz', '', filename)
    fields = filename.split('.')
    ok = True
    if (len(fields) != 11 and len(fields) != 19):
        print(str(len(fields)) + ' fields\t' + uri)
        return None
    try:
        freeze = int(fields[3])
        n_cases = int(fields[7])
        n_controls = int(fields[8])
        date = int(fields[10])
    except ValueError:
        print('Number error\t{}'.format(uri))
        ok = False
    pop = fields[6].upper()
    if pop not in POP:
        print('Unexpected pop {}\t{}'.format(pop, uri))
        ok = False
    ana = None
    for a in ANA:
        if a in fields[2]:
            ana = a
            break
    if ana is None:
        print('Unexpected ANA {}\t{}'.format(fields[2], uri))
        ok = False
    age = fields[4].upper()
    if age not in AGE:
        print('Unexpected age {}\t{}'.format(age, uri))
        ok = False
    s = fields[5].upper()
    if s.startswith('F'):
        sex = 'FEMALE'
    elif s.startswith('M'):
        sex = 'MALE'
    elif s == 'ALL':
        sex = 'ALL'
    else:
        print('Unexpected sex {}\t{}'.format(fields[5], uri))
        ok = False
    if not ok:
        return None
    return {'cohort': fields[0], 'file': uri, 'ana': ana, 'age': age, 'sex': sex, 'pop': pop, 'n_cases': n_cases, 'n_controls': n_controls}

def check_files(filename):
    with open(filename, 'rt') as f:
        parsed = [parse_file(line.strip()) for line in f if not line.startswith('#')]
    for p in parsed:
        if p is None:
            raise(Exception('File check failed'))
    return parsed

def write_munge_input(studies):
    allowed_pops = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas', 'all']
    replace_pops = [
        ('non-eur', 'all'),
        ('arab', 'all'),
        ('eur', 'nfe'),
        ('csa', 'sas'),
        ('his', 'amr')
    ]
    out = open('munge_input.txt', 'wt')
    for study in studies:
        pop_gn = study['pop'].lower()
        for pop in replace_pops:
            if pop_gn == pop[0]:
                pop_gn = pop[1]
                break
        if pop_gn not in allowed_pops:
            raise(Exception('Unexpected population {} in {}'.format(pop_gn, study['file'])))
        out.write('\t'.join([study['file'], pop_gn, str(study['n_cases'] + study['n_controls'])]) + '\n')
    out.close()

def create_meta_conf(studies):
    conf = dd(list)
    for study in studies:
        cohort = {
            'name': study['cohort'] + '_' + study['pop'],
            'file': study['file'].replace('gs://', '/cromwell_root/'),
            'n_cases': study['n_cases'],
            'n_controls': study['n_controls'],
            'chr':'#CHR',
            'pos':'POS',
            'ref':'Allele1',
            'alt':'Allele2',
            'effect':'BETA',
            'effect_type':'beta',
            'pval':'p.value',
            'se':'SE',
            'extra_cols':['AF_Allele2', 'AF_fc', 'imputationInfo', 'N']
        }
        if study['sex'] != 'ALL' and study['age'] != 'ALL':
            #print('Both age and sex stratified\t{}'.format(study['file']))
            continue
        if study['sex'] == 'ALL':
            analysis = study['ana'] + '_' + study['age']
        elif study['age'] == 'ALL':
            analysis = study['ana'] + '_' + study['sex']
        else:
            raise Exception('Unexpected study {}'.format(study))
        conf[analysis].append(cohort)

    for ana in [('A1', 'A2'), ('B1', 'B2'), ('C1', 'C2'), ('A1', 'B1'), ('A2', 'B2'), ('B2', 'C2')]:
        for type in ['ALL', 'MALE', 'FEMALE', 'GT_60', 'LE_60']:
            ana1 = ana[0] + '_' + type
            ana2 = ana[1] + '_' + type
            if ana1 not in conf or ana2 not in conf:
                continue
            for cohort in conf[ana1]:
                if cohort['name'] not in [c['name'] for c in conf[ana2]]:
                    conf[ana2].append(cohort)
    # remove B1 from C2 (as per above it's possible to go from B1 to B2 and then from B2 to C2 - but in C2, all B1 individuals should be cases)
    # nice code
    remove = []
    for ana1 in conf:
        if ana1.startswith('B1'):
            for ana2 in conf:
                if ana2.startswith('C2'):
                    for c in conf[ana2]:
                        if c in conf[ana1]:
                            remove.append((ana2, c['name']))
    for r in remove:
        conf[r[0]] = [c for c in conf[r[0]] if c['name'] != r[1]]

    #remove double HOSTAGE
    conf['C2_ALL'] = [c for c in conf['C2_ALL'] if not (c['name'].startswith('Italy_HOSTAGE') or c['name'].startswith('Spain_HOSTAGE'))]
    conf['B2_ALL'] = [c for c in conf['B2_ALL'] if not (c['name'].startswith('Italy_HOSTAGE') or c['name'].startswith('Spain_HOSTAGE'))]

    out = open('analysis_files.txt', 'wt')
    out_pheno = open('analysis_phenos.txt', 'wt')
    if not os.path.exists('json'):
        os.makedirs('json')
    if not os.path.exists('json/leave'):
        os.makedirs('json/leave')
    for analysis in od(sorted(conf.items())):
        cohorts = set([cohort['name'].split('_')[0] for cohort in conf[analysis]])
        for c in cohorts:
            leave = [cohort for cohort in conf[analysis] if not cohort['name'].startswith(c)]
            with open('json/leave/{}_leave_{}.json'.format(analysis, c), 'wt') as f:
                json.dump({'meta': leave}, f, indent=4)
            if analysis == 'B2_ALL' or analysis == 'C2_ALL':
                if analysis == 'B2_ALL':
                    n=10
                if analysis == 'C2_ALL':
                    n=12
                out_pheno.write('{}_leave_{}\t{}\n'.format(analysis, c, n))
    quit()
    for analysis in od(sorted(conf.items())):
        print(analysis)
        with open('json/{}.json'.format(analysis), 'wt') as f:
            json.dump({'meta': conf[analysis]}, f, indent=4)
        leave_ukbb = [cohort for cohort in conf[analysis] if not cohort['name'].lower().startswith('ukbb')]
        with open('json/{}_leave_UKBB.json'.format(analysis), 'wt') as f:
            json.dump({'meta': leave_ukbb}, f, indent=4)
        leave_genomicc = [cohort for cohort in conf[analysis] if not cohort['name'].lower().startswith('genomicc')]
        with open('json/{}_leave_genomicc.json'.format(analysis), 'wt') as f:
            json.dump({'meta': leave_genomicc}, f, indent=4)
        leave_23andme = [cohort for cohort in conf[analysis] if not cohort['name'].lower().startswith('23andme')]
        with open('json/{}_leave_23andme.json'.format(analysis), 'wt') as f:
            json.dump({'meta': leave_23andme}, f, indent=4)
        leave_prs = [cohort for cohort in conf[analysis] if not (
            cohort['name'].lower().startswith('hostage') or
            cohort['name'].lower().startswith('belcovid') or
            cohort['name'].lower().startswith('bosco') or
            cohort['name'].lower().startswith('gencovid') or
            cohort['name'].lower().startswith('bracovid') or
            cohort['name'].lower().startswith('swecovid') or
            cohort['name'].lower().startswith('bqc19') or
            cohort['name'].lower().startswith('spgrx')
            )]
        with open('json/{}_leave_prs.json'.format(analysis), 'wt') as f:
            json.dump({'meta': leave_prs}, f, indent=4)
        admixed = [cohort for cohort in conf[analysis] if cohort['name'].endswith('_AFR') or cohort['name'].endswith('_HIS') or cohort['name'].endswith('AMR')]
        with open('json/{}_admixed.json'.format(analysis), 'wt') as f:
            json.dump({'meta': admixed}, f, indent=4)
        eur = [cohort for cohort in conf[analysis] if cohort['name'].endswith('_EUR') or cohort['name'].endswith('_FIN')]
        with open('json/{}_eur.json'.format(analysis), 'wt') as f:
            json.dump({'meta': eur}, f, indent=4)
        eur_leave23andme = [cohort for cohort in conf[analysis] if (
            cohort['name'].endswith('_EUR') or
            cohort['name'].endswith('_FIN')
            ) and not (
            cohort['name'].lower().startswith('23andme')
            )]
        with open('json/{}_eur_leave_23andme.json'.format(analysis), 'wt') as f:
            json.dump({'meta': eur_leave23andme}, f, indent=4)
        eur_leaveukbb = [cohort for cohort in conf[analysis] if (
            cohort['name'].endswith('_EUR') or
            cohort['name'].endswith('_FIN')
            ) and not (
            cohort['name'].lower().startswith('ukbb')
            )]
        with open('json/{}_eur_leave_ukbb.json'.format(analysis), 'wt') as f:
            json.dump({'meta': eur_leaveukbb}, f, indent=4)
        eur_leaveukbb23andme = [cohort for cohort in conf[analysis] if (
            cohort['name'].endswith('_EUR') or
            cohort['name'].endswith('_FIN')
            ) and not (
            cohort['name'].lower().startswith('ukbb') or
            cohort['name'].lower().startswith('23andme')
            )]
        with open('json/{}_eur_leave_ukbb_23andme.json'.format(analysis), 'wt') as f:
            json.dump({'meta': eur_leaveukbb23andme}, f, indent=4)
        for cohort in conf[analysis]:
            print('{}\t{}\t{}'.format(cohort['name'], cohort['n_cases'], cohort['n_controls']))
            out.write(cohort['file'].replace('/cromwell_root/', 'gs://') + '\t')
        out.write('\n')

        print()
    out.close()

if __name__ == '__main__':
    studies = check_files(sys.argv[1])
    write_munge_input(studies)
    create_meta_conf(studies)
