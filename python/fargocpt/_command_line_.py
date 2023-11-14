import argparse

def main():

    parser = argparse.ArgumentParser(description='FargoCPT: A Python package for handling FargoCPT simulations.')
    parser.add_argument('command', type=str, help='Command to execute.', choices=["run", "data"])
    opts, remainder = parser.parse_known_args()

    print(opts, remainder)

    if opts.command == "run":
        from .run import main as cmd_main
        cmd_main(remainder)
    elif opts.command == "data":
        from .data import main as cmd_main
        cmd_main(remainder)

if __name__=="__main__":
    main()