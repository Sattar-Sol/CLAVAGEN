"""Where CLAVAGEN keeps its settings, reference library and database cache.

The folder is ~/.clavagen. Earlier releases (ACMG-VarCall) used ~/.acmg_varcall: on first start
its contents are copied over, so saved settings, the reference library and the cache are kept.
"""
import os
import shutil

OLD_HOME = os.path.join(os.path.expanduser("~"), ".acmg_varcall")


def app_home():
    home = os.environ.get("CLAVAGEN_HOME") or os.path.join(os.path.expanduser("~"), ".clavagen")
    if not os.path.isdir(home):
        if os.path.isdir(OLD_HOME) and not os.environ.get("CLAVAGEN_HOME"):
            try:
                shutil.copytree(OLD_HOME, home, ignore=shutil.ignore_patterns("window_profile"))
            except OSError:
                os.makedirs(home, exist_ok=True)
        else:
            os.makedirs(home, exist_ok=True)
    return home


def app_path(*parts):
    return os.path.join(app_home(), *parts)
