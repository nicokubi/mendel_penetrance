import argparse
import pexpect
import sys

def run_mendel(mendel='single', cwd=None):
    proc = pexpect.spawn([mendel], cwd=cwd)
    # proc.logfile = sys.stdout.buffer
    # DO YOU WISH TO CHANGE TO BATCH MODE? [YES/NO]
    proc.expect(':')
    proc.sendline('no')
    # CHOOSE AN ITEM [1,...,21]:
    proc.expect(':')
    proc.sendline('21')
    # ANOTHER PROBLEM [YES/NO]:
    proc.expect(':')
    proc.sendline('no')
    proc.interact()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mendel', default='../mendel/single', help='full path to single command')
    parser.add_argument('--cwd', help='Current working directory')
    args = parser.parse_args()

    run_mendel(mendel=args.mendel, cwd=args.cwd)
