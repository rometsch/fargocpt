import argparse
import sys

subcommands = ["run", "data", "config"]

def main(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description='FargoCPT: A Python package for handling FargoCPT simulations.')
    parser.add_argument('command', type=str, help='Command to execute.', choices=subcommands)

    # handle -h and --help flags for subcommands
    # pass on to subcommand if -h/--help comes after subcommand
    pass_on_help = False
    if "-h" in args or "--help" in args:
        ind = args.index("-h") if "-h" in args else args.index("--help")
        first_sc_ind = 1000
        for sc in subcommands:
            if sc in args:
                first_sc_ind = min(args.index(sc), first_sc_ind)
        if ind > first_sc_ind:
            # pass on help to subcommand
            pass_on_help = True
        # remove help from arguments
        args.pop(ind)

    opts, remainder = parser.parse_known_args(args)

    if pass_on_help:
        remainder.append("--help")

    if opts.command == "run":
        from .run import main as cmd_main
        cmd_main(remainder)
    elif opts.command == "data":
        from .data import main as cmd_main
        cmd_main(remainder)
    elif opts.command == "config":
        from .config import main as cmd_main
        cmd_main(remainder)

if __name__=="__main__":
    main()