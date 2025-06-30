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
