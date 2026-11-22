"""Compatibility alias: ACMG-VarCall is now CLAVAGEN.

``python3 -m acmg_varcall ...`` and ``import acmg_varcall.<module>`` keep working and use the
``clavagen`` package. New code should import ``clavagen``.
"""
import importlib
import os
import importlib.abc
import importlib.util
import sys

import clavagen
from clavagen import __version__  # noqa: F401


class _AliasFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def find_spec(self, name, path=None, target=None):
        if name.startswith("acmg_varcall.") and name != "acmg_varcall.__main__":
            return importlib.util.spec_from_loader(name, self)
        return None

    def create_module(self, spec):
        return importlib.import_module("clavagen." + spec.name.split(".", 1)[1])

    def exec_module(self, module):
        pass


sys.meta_path.insert(0, _AliasFinder())
__path__ = [os.path.dirname(os.path.abspath(__file__))]  # only __main__.py lives here; the rest is aliased
