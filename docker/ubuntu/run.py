import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument("--arg-1", help="an arg")
parser.add_argument("--arg-2", help="another arg")
args = parser.parse_args()
#print(subprocess.check_output(["echo", "Testing output; args:", args.arg_1, args.arg_2]))
print("Testing output; args:", args.arg_1, args.arg_2)

