import argparse

parser = argparse.ArgumentParser()

# An int is an explicit number of arguments to accept.
parser.add_argument('--nargs', nargs='+')

# To make the input integers
parser.add_argument('-a','--nargs-int-type', default=[0], nargs='+', type=int)

# An alternate way to accept multiple inputs, but you must
# provide the flag once per input. Of course, you can use
# type=int here if you want.
parser.add_argument('--append-action', action='append')

print(parser.parse_args().nargs_int_type)
# To show the results of the given option to screen.
for ARG, value in parser.parse_args()._get_kwargs():
    if value is not None:
        print(ARG,value)