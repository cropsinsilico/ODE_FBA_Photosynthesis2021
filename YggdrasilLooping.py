from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args = parser.parse_args(["yggdrasil_ODE_FBA_testing.yaml"])
runner.run(args.yamlfile, ygg_debug_prefix=prog)
