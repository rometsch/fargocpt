import yaml
from tabulate import tabulate

def convert_to_md(datafile, outfile):
    # Load the YAML data
    with open(datafile, 'r') as infile:
        parsed_data = yaml.safe_load(infile)

    # Convert the parsed data into a list of lists for the tables
    table_data = [["Parameter", "Choices", "Default", "Type", "Unit Support", "Description"]]
    deprecated_data = [["Parameter", "Newname", "Hint", "Choices", "Default", "Type", "Unit Support", "Description"]]
    for key, value in parsed_data.items():
        deprecated = value.get('deprecated', '')
        if not deprecated:
            row = [
                key,
                value.get('choices', ''),
                value.get('default', ''),
                value.get('type', ''),
                value.get('unitsupport', ''),
                value.get('description', '')
            ]
            table_data.append(row)
        else:
            row = [
                key,
                value.get('newname', ''),
                value.get('hint', ''),
                value.get('choices', ''),
                value.get('default', ''),
                value.get('type', ''),
                value.get('unitsupport', ''),
                value.get('description', '')
            ]
            deprecated_data.append(row)

    table_header = table_data[0]
    deprecated_header = deprecated_data[0]
    # Sort the table data by parameter name
    table_data = sorted(table_data[1:], key=lambda row: row[0])
    deprecated_data = sorted(deprecated_data[1:], key=lambda row: row[0])

    # Add the headers back to the start of the table data
    table_data.insert(0, table_header)
    deprecated_data.insert(0, deprecated_header)

    # Generate the markdown tables
    markdown_table = tabulate(table_data, headers="firstrow", tablefmt="pipe")
    deprecated_table = tabulate(deprecated_data, headers="firstrow", tablefmt="pipe")


    with open(outfile, 'w') as outfile:
        outfile.write("# Parameters\n")
        outfile.write("Following is a list of parameters that can be used to configure physics and numerics of a simulation and the behavior of the program. \n\n")
        outfile.write("## Active Parameters\n")
        outfile.write(" For numerical parameters, the range of choices are abreviated as follows:\n")
        outfile.write("- '-' : the range includes negative numbers.\n")
        outfile.write("- '0' : the range includes zero.\n")
        outfile.write("- '+' : the range includes positive numbers.\n\n")
        outfile.write(markdown_table)
        outfile.write("\n\n## Deprecated Parameters\n")
        outfile.write(deprecated_table)


if __name__ == "__main__":
    # parse optional cli args for input and output files
    import argparse
    parser = argparse.ArgumentParser(description='Convert a YAML file to a markdown table.')
    parser.add_argument('-i','--input', default="parameters.yml", help='The YAML file to convert.')
    parser.add_argument('-o','--output', default="parameters.md", help='The markdown file to write to.')
    args = parser.parse_args()


    convert_to_md(args.input, args.output)