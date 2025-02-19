import os
import pandas as pd
import glob
import yaml
import argparse

# Define a custom list class to force inline (flow style) representation.
class FlowSeq(list):
    pass

# Register a YAML representer for FlowSeq so that the list is output inline.
def flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(FlowSeq, flow_seq_representer)

# Extracts the full basename
def extract_basename(file_path):
    """
    Returns basename:
      'BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP'
    from a path:
      '/.../BCA001_lib_13077AAF_CAGATCAC-ATGTGAAG_DSP/Solo.out/Gene/UMIperCellSorted.txt'
    """
    base = os.path.basename(os.path.dirname(file_path))
    config = os.path.basename(os.path.dirname(os.path.dirname(file_path))).rsplit('_')[1]

    return base + '_' + config


# Loop through all subdirectories in the saturation directory and create dictionary
def get_data(all_files):
    data_dict = {}
    
    for file in all_files:
        sample_name = extract_basename(file)
        df = pd.read_csv(file, sep='\t')
        data_dict[sample_name] = df[['ninput', 'sat']].to_dict(orient='list')
    print(data_dict)
    return data_dict



def main(baseDir):
    # Retrieve paths for all saturation_output.tsv files found in base directory
    # example: baseDir = "/no_backup/asebe/bvanwaardenburg/data/250115_ParseBio_Nvec_Tcas_Pliv_Cele/Nvec_BCA009_BCA010"
    pattern = "**/saturation_output.tsv"
    all_files = glob.glob(os.path.join(baseDir, pattern), recursive=True)

    # Saves ninput and saturation values for each file
    data = get_data(all_files)

    # Transform the data: wrap the lists in FlowSeq to force inline YAML representation.
    transformed_data = {"data": []}
    for key, values in data.items():
        transformed_data["data"].append({
            "label": f"bca_{key}",
            "x": FlowSeq(values.get("ninput", [])),
            "y": FlowSeq(values.get("sat", []))
        })


    # Create the YAML content combining header details and data
    content = {
        "id": "saturation_lineplot",
        "section_name": "Saturation Analysis",
        "description": "This plot shows the saturation analysis results.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": "saturation_lineplot",
            "title": "Saturation Analysis",
            "xlab": "Number of Input Reads",
            "ylab": "Saturation",
            "ymin": 0,
            "ymax": 1,
            "xLog": True
        },
        "data": transformed_data
    }

    # Convert to YAML format and save
    with open("%s/saturation_comb_mqc.yaml" %(baseDir), "w") as f:
        yaml.dump(content, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combines saturation output data into YAML file.")
    parser.add_argument("baseDir", help="Base directory path.")
    args = parser.parse_args()

    main(args.baseDir)