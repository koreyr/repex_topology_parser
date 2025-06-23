#!/bin/env python3
import argparse
parser = argparse.ArgumentParser(description="Flip a switch by setting a flag")
parser.add_argument('-w', action='store_true')

args = parser.parse_args()
print(args.w)
