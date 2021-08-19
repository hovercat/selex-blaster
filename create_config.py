#!/usr/bin/env python
import os


######################################################
# Argument parsing functions
######################################################
def parse_string(s, default=None):
    s = s.strip()
    if s is not None and s != "":
        return s
    else:
        return default


def parse_list(s, delimiter=' ', default=None):
    s = s.strip()
    return s.split(delimiter)


def parse_boolean(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s) == 1
    else:
        return default


def parse_float(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return float(s)
    else:
        return default


def parse_int(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s)
    else:
        return default


def make_sane_string(s):
    return ''.join([char for char in s.strip() if char.isalpha() or char.isdigit() or char == ' ']).rstrip()


######################################################
# Starting config creation
######################################################
arguments = dict()

print("Creating a config file for selex-blaster workflow.")
print(
    "Please answer the following questions. If you don't know an answer, you can press enter without an answer to opt for the default value.")

# Experiment Name
while not arguments.get("experiment"):
    arguments["experiment"] = parse_string(input("What's your SELEX experiment's unique identifying name?\n"))
    if not arguments["experiment"]:
        print("There is no default value for the experiment's name.")

# Input Directory
arguments["input_dir"] = parse_string(input("In which directory are the preprocessed FASTA-files stored? (default: "
                                            "./input_dir)\n"), "./input_dir")

# Output Directory
arguments["output_dir"] = parse_string(
    input("In which directory to write output files to? (default: './output_{}/)\n".format(make_sane_string(arguments["experiment"]))),
    './output_{}'.format(make_sane_string(arguments["experiment"]))
)

# Check input directory for FASTA-Files
if os.path.exists(arguments["input_dir"]):
    fasta_files = [f for f in os.listdir(arguments["input_dir"]) if
                   os.path.isfile(os.path.join(arguments["input_dir"], f))
                   and (f.lower().endswith(".fasta") or f.lower().endswith(".fa"))]
    fasta_files.sort()

    print("Files found in the input directory '{}'".format(arguments["input_dir"]))
    print('\n'.join(fasta_files))
else:
    print("Input directory '{}' does not exist.".format(arguments["input_dir"]))

# Ask for round names
while not arguments.get("round_order"):
    arguments["round_order"] = parse_list(
        input(
            "Please provide the list of SELEX rounds in sequential order, separated by a space. (e.g. R0 R2 R4 R6)\n"),
        delimiter=" ")
    if not arguments.get("round_order"):
        print("Please specify the round order of your FASTA-files!")
while not arguments.get("library"):
    arguments["library"] = parse_string(input("Which SELEX round is the library?"))


# Ask for FASTA-Pattern
while not arguments.get("fasta_pattern"):
    arguments["fasta_pattern"] = parse_string(input(
        "Please provide a FASTA-search pattern, including a wild-card. For example: \nIf you had R0.fasta your search pattern would be: \'*1.fasta\'\n"))
    if not arguments["fasta_pattern"]:
        print("Please specify a FASTA-search pattern!")

while not arguments.get("primer5"):
    arguments["primer5"] = parse_string(input("Please provide the forward primer (5'-primer).\n"))
    if not arguments["primer5"]:
        print("Please specify a forward primer!")

while not arguments.get("primer3"):
    arguments["primer3"] = parse_string(input("Please provide the reverse primer (3'-primer).\n"))
    if not arguments["primer3"]:
        print("Please specify a reverse primer!")

while not arguments.get("folding_temp"):
    arguments["folding_temp"] = parse_int(input("What's the temperature used during your experiment? Provide a whole number.\n"))
while not arguments.get("cpus"):
    arguments["cpus"] = parse_int(input("What's your max number of available CPUs?\n"))
while not arguments.get("inflation"):
    arguments["inflation"] = parse_float(input("Provide a inflation value for MCL. default=1.4\n"), 1.4)
while not arguments.get("subsample_library"):
    arguments["subsample_library"] = parse_float(input("How many sequences should be resampled from the library for motif detection? default: 10000\n"), 10000)
while not arguments.get("RNAfold_mathews2004_dna"):
    arguments["RNAfold_mathews2004_dna"] = parse_string(input("Please provide the location where Mathews DNA 2004 model is.\n"))


# Finished reading config parameters.
print("")
config_file_name = "./{}.config".format(make_sane_string(arguments["experiment"]))
with open(config_file_name, "w") as config_file:
    config_file.write("params {\n")
    for key, value in arguments.items():
        config_file.write("\t")
        if isinstance(value, list):
            config_file.write('{key} = ["'.format(key=key))
            config_file.write('","'.join(value))
            config_file.write('"]\n')
        elif isinstance(value, int):
            config_file.write('{key} = {value}\n'.format(key=key, value=value))
        elif isinstance(value, str):
            config_file.write('{key} = "{value}"\n'.format(key=key, value=value))
        elif isinstance(value, bool):
            if value:
                config_file.write('{key} = true\n')
            else:
                config_file.write('{key} = false\n')
    config_file.write("}\n")

print("Done.")
print("Your config file was saved to '{}'.".format(config_file_name))
