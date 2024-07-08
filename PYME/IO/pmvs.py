from PYME.warnings import warn

def chkpath_relative(pmvsname,filename):
    from pathlib import Path
    pf = Path(filename)
    if pf.is_absolute():
        return filename
    parentdir = Path(pmvsname).parent
    if parentdir == Path('.'):
        return filename
    return str(parentdir / filename)

def check_entries(pmvs,required=[],optional=[],deprecated={}):
    required_flat = []
    for entry in required:
        if isinstance(entry,(list,tuple)):
            missing = True
            for subentry in entry:
                if subentry in pmvs:
                    missing = False
                required_flat.append(subentry)
            if missing:
                raise RuntimeError("fatal: none of required '%s' entry alternatives in PMVS file" % entry)
        else:
            if entry not in pmvs:
                raise RuntimeError("fatal: no required '%s' entry in PMVS file" % entry)
            required_flat.append(entry)
    for entry in pmvs:
        if entry not in required_flat + optional:
            raise RuntimeError("fatal: unknown '%s' entry in PMVS file" % entry)
    for entry in deprecated.keys():
        if entry in pmvs:
            warn("deprecated key '%s' in pmvs file, replace with current version '%s'" % (entry,deprecated[entry]))

def load_pmvs(filename,translate_paths=True):
    import yaml
    with open(filename, 'r') as file:
        pmvs_args = yaml.safe_load(file)
    if pmvs_args.get('pmvs_version') != 'v1.0':
        raise RuntimeError("file %s is not a PMVS version 1.0 file" % filename)
    check_entries(pmvs_args,
                  required=[['localizations','mainfile'],'pmvs_version'],
                  optional=['load','imageds','recipe','comment'],
                  deprecated=dict(localizations='mainfile',imageds='load'))
    if 'imageds' in pmvs_args:
        # use imageds only for backwards compatibility
        if 'load' not in pmvs_args:
            pmvs_args['load'] = {}
        for key in pmvs_args['imageds']:
            pmvs_args['load'][key] = pmvs_args['imageds'][key]
    if 'localizations' in pmvs_args:
        pmvs_args['mainfile'] = pmvs_args['localizations']
    if not translate_paths:
        return pmvs_args
    
    # otherwise make sure relative paths are suitably translated
    for entry in ['recipe', 'mainfile']:
        if entry in pmvs_args:
            pmvs_args[entry] = chkpath_relative(filename,pmvs_args[entry])
    if 'load' in pmvs_args:
        for key in pmvs_args['load']:
            pmvs_args['load'][key] = chkpath_relative(filename,pmvs_args['load'][key])
    return pmvs_args

def load_and_parse_pmvs(args):
    pmvs_args = load_pmvs(args.file)
    if len(args.load) > 0: # check if this restriction is really needed
        raise RuntimeError("loading additional files from the command line not allowed when using pmvs file")
    if pmvs_args.get('recipe',None) is not None:
        if args.recipe is not None: # check if restriction sensible
            raise RuntimeError("recipe in pmvs file but recipe also supplied on command line - conflict")
        args.recipe = pmvs_args['recipe']
    for name in pmvs_args.get('load',{}):
        args.load.append((name,pmvs_args['load'][name]))
    args.file = pmvs_args['mainfile']
    return args
