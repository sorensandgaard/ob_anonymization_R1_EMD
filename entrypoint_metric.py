import argparse
import os

def concatenate_input_content(input_files):
    concatenated_content = ""  # Initialize an empty string to hold the concatenated content
    
    # Iterate over each input file
    for input_file in input_files:
        # Open each input file in read mode and read its content
        with open(input_file, 'r') as file:
            # Read the content of the input file and append it to the concatenated_content string
            concatenated_content += file.read()
            # Optionally, you can add a newline between the content of each file
            concatenated_content += '\n'

    return concatenated_content


def run_metric(output_dir, name, input_files):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    content = concatenate_input_content(input_files)

    metric_results_file = os.path.join(output_dir, f'{name}.results.txt')
    content += f"4. Running metric using parameters into {metric_results_file}"

    with open(metric_results_file, 'w') as file:
        file.write(content)


def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Run metic on files.')

    # Add arguments
    parser.add_argument('--output_dir', type=str, help='output directory where metic will store results.')
    parser.add_argument('--name', type=str, help='name of the dataset')
    parser.add_argument('--methods.mapping', type=str, help='input file #1.')
    parser.add_argument('--data.meta', type=str, help='input file #2.')
    parser.add_argument('--data.data_specific_params', type=str, help='input file #3.')

    # Parse arguments
    args = parser.parse_args()

    methods_mapping_input = getattr(args, 'methods.mapping')
    data_meta_input = getattr(args, 'data.meta')
    data_params_input = getattr(args, 'data.data_specific_params')
    input_files = [methods_mapping_input, data_meta_input, data_params_input]

    run_metric(args.output_dir, args.name, input_files)


if __name__ == "__main__":
    main()
