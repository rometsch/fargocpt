import yaml
from tabulate import tabulate

# Load the YAML data
with open("parameters.yml", 'r') as infile:
    parsed_data = yaml.safe_load(infile)

# Convert the parsed data into a list of lists for the tables
table_data = [["Parameter", "Choices", "Default", "Description", "Type", "Unit Support", "Deprecated", "Hint"]]
deprecated_data = [["Parameter", "Choices", "Default", "Description", "Type", "Unit Support", "Deprecated", "Hint"]]
for key, value in parsed_data.items():
    row = [
        key,
        value.get('choices', ''),
        value.get('default', ''),
        value.get('description', ''),
        value.get('type', ''),
        value.get('unitsupport', ''),
        value.get('deprecated', ''),
        value.get('hint', '')
    ]
    if value.get('deprecated') == True:
        deprecated_data.append(row)
    else:
        table_data.append(row)

# Sort the table data by parameter name
table_data = sorted(table_data[1:], key=lambda row: row[0])
deprecated_data = sorted(deprecated_data[1:], key=lambda row: row[0])

# Add the headers back to the start of the table data
table_data.insert(0, ["Parameter", "Choices", "Default", "Description", "Type", "Unit Support", "Deprecated", "Hint"])
deprecated_data.insert(0, ["Parameter", "Choices", "Default", "Description", "Type", "Unit Support", "Deprecated", "Hint"])

# Generate the markdown tables
markdown_table = tabulate(table_data, headers="firstrow", tablefmt="pipe")
deprecated_table = tabulate(deprecated_data, headers="firstrow", tablefmt="pipe")


with open("parameters.md", 'w') as outfile:
    outfile.write("# Active Parameters\n")
    outfile.write(" For numerical parameters, the range of choices are abreviated as follows:\n")
    outfile.write("- '-' : the range includes negative numbers.\n")
    outfile.write("- '0' : the range includes zero.\n")
    outfile.write("- '+' : the range includes positive numbers.\n\n")
    outfile.write(markdown_table)
    outfile.write("\n\n# Deprecated Parameters\n")
    outfile.write(deprecated_table)