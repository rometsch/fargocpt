""" Config structure for fargocpt python wrapper. """
import os
import sys
import yaml

program_name = "fargocpt"
config_version = "1.0"
information_types = ["exe_path"]

config_dir = os.path.join(os.path.expanduser("~"), f".config/{program_name}")

def main(args=sys.argv[1:]):
    opts = parse_opts(args)
    c = Config()

    if opts.subparser_name in ["show"] or opts.subparser_name is None:
        opts.func(c)
    elif opts.subparser_name in ["get"]:
        opts.func(c, opts.key)
    else:
        opts.func(c, opts.key, opts.value)


def parse_opts(args=sys.argv[1:]):
    import argparse
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=Config.print)

    subparsers = parser.add_subparsers(dest='subparser_name')
    parser_set = subparsers.add_parser('set', help='Set a config item.')
    parser_set.add_argument("key",
                            choices=information_types,
                            help="What to set.")
    parser_set.add_argument("value")
    parser_set.set_defaults(func=Config.set)

    parser_remove = subparsers.add_parser('remove',
                                          help='Remove a config item.')
    parser_remove.add_argument("key",
                               choices=information_types,
                               help="What to set.")
    parser_remove.add_argument("value")
    parser_remove.set_defaults(func=Config.remove)

    parser_show = subparsers.add_parser('show', help='Show the config.')
    parser_show.set_defaults(func=Config.print)
    parser_get = subparsers.add_parser(
        'get', help='Return the value of a root level config item.')
    parser_get.add_argument("key", help="What to get.")
    parser_get.set_defaults(func=Config.print_value)
    args = parser.parse_args(args)
    return args


def expand_path(path):
    abspath = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(abspath):
        raise FileNotFoundError("No such directory: {}".format(path))
    return abspath

def check_information_type(info_type):
    if not any((info_type == t for t in information_types)):
        raise AttributeError(
            "Information type {} not supported".format(info_type))


class Config:
    def __init__(self):
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
        self.config_file = os.path.join(config_dir, "config.yml")
        self.load()

    def set(self, what, val):
        check_information_type(what)
        # use absolute path if val is a path
        if what == "exe_path":
            val = expand_path(val)
        self.data[what] = val
        self.save()

    def remove(self, what):
        check_information_type(what)
        list_name = what + "_list"
        try:
            del self.data[what]
            self.save()
        except KeyError:
            print("No config for type", what)
            pass

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, val):
        self.data[key] = val

    def save(self):
        self.data["type"] = f"{program_name} config"
        self.data["version"] = config_version
        with open(self.config_file, "w") as outfile:
            yaml.dump(self.data, outfile, default_flow_style=False)

    def load(self):
        try:
            with open(self.config_file, "r") as infile:
                self.data = yaml.safe_load(infile)
        except FileNotFoundError:
            self.data = {}
            self.data["exe_path"] = None

    def print(self):
        try:
            from pprint import pprint as print
        except ImportError:
            pass
        print(self.data)

    def print_value(self, key):
        try:
            print(self[key])
        except KeyError:
            print("Error: No config value found for key '{}'".format(key))
            sys.exit(1)
