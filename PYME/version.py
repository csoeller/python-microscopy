
# PYME uses date based versions (yy.mm.dd)
# NOTE: This file is automatically generated as part of the CI process. DO NOT EDIT. Changes should be made to the temlate in misc/update_version.py instead.
_release_version = '25.03.27'
_release_changeset = '864f1f84cce3'

version = _release_version
changeset = _release_changeset # Git changeset ID

# if we are a development install, modify our version number to indicate this
# note this duplicates the logic in PYME.misc.check_for_updates, but is reproduced verbatim here
# so that the version.py works even if PYME is not already installed (ie when you are running python setup.py develop)
import os
pyme_parent_dir = os.path.dirname(os.path.dirname(__file__))
    
if os.path.exists(os.path.join(pyme_parent_dir, '.git')):
    print('Detected a .git folder, assuming a development install')
    dev_install = True
else:
    dev_install = False

if dev_install:
    # development install, flag this as a postrelease (means it should be given higher priority than an conda package of the same release version)
    version = version + '.post0.dev'

    try:
        import subprocess

        p = subprocess.Popen('git describe --abbrev=12 --always --dirty=+', shell=True, stdout=subprocess.PIPE, encoding='utf8')
        changeset = p.stdout.readline().strip()
    except:
        print('Development install detected, but could not get git hash, is git installed?')


_detailed_version = None
def detailed_version():
    '''
    Version to display in about dialogs, error dialogs, etc .... includes install type and commit hash (if modified from release commit). 
    
    Is NOT pep-0440 compliant, as it requires human interpretation of, e.g. commit hash ordering. 

    Example full_version strings:

    21.10.01[conda] - a conda package based install from an official release (also executable installers)
    21.10.01[pip] - a pip install from an official release
    21.10.01.post0.dev[git] - a development install of the exact release version
    21.10.01.post0.dev[git]f94cc30be308 - a development install which has been modified since the release (note, does not distinguish between remote commits to master and local commits)
    21.10.01.post0.dev[git]f94cc30be308+ - a development install with uncommitted local changes

    '''
    global _detailed_version
    if _detailed_version is None:
        from PYME.misc.check_for_updates import guess_install_type

        fv = version + '[' + guess_install_type() + ']'
        if changeset !=_release_changeset:
            # code has been modified since last release, append commit hash
            fv += changeset

        _detailed_version = fv

    return _detailed_version
