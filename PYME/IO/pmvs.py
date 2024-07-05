def chkpath_relative(pmvsname,filename):
    from pathlib import Path
    pf = Path(filename)
    if pf.is_absolute():
        return filename
    parentdir = Path(pmvsname).parent
    if parentdir == Path('.'):
        return filename
    return str(parentdir / filename)

def check_entries(pmvs,required=[],optional=[]):
    for entry in required:
        if entry not in pmvs:
            raise RuntimeError("fatal: no required '%s' entry in PMVS file" % entry)
    for entry in pmvs.keys():
        if entry not in required + optional:
            raise RuntimeError("fatal: unknown '%s' entry in PMVS file" % entry)

def load_pmvs(filename,translate_paths=True):
    import yaml
    with open(filename, 'r') as file:
        pmvs_args = yaml.safe_load(file)
    if pmvs_args.get('pmvs_version') != 'v1.0':
        raise RuntimeError("file %s is not a PMVS version 1.0 file" % filename)
    check_entries(pmvs_args,
                  required=['localizations','pmvs_version'],
                  optional=['imageds','recipe','comment'])
    
    if not translate_paths:
        return pmvs_args
    
    # otherwise make sure relative paths are suitably translated
    for entry in ['localizations', 'recipe']:
        if entry in pmvs_args:
            pmvs_args[entry] = chkpath_relative(filename,pmvs_args[entry])
    for key in pmvs_args['imageds']:
        pmvs_args['imageds'][key] = chkpath_relative(filename,pmvs_args['imageds'][key])
    return pmvs_args
