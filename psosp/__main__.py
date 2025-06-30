"""
This module serves as the command-line entry point for the psosp package.

It is created to avoid the name collision issue where the executable script `psosp`
conflicts with the package name `psosp`, causing a `ModuleNotFoundError` during
execution in certain environments like conda-build.

This setup allows the tool to be run via `psosp` command or `python -m psosp`.
"""
import sys
from . import main as psosp_main

def main():
    """
    Wrapper for the main function.
    This is the function that the `entry_points` in setup.py will call.
    """
    sys.exit(psosp_main())

if __name__ == '__main__':
    main() 