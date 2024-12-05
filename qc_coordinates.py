import pandas as pd

input_file = snakemake.input.data[0]
coordinates_span = snakemake.params.coordinates_span
output_file = snakemake.output[0]
    
data = pd.read_csv(input_file, dtype={'device_id': str, 'subject_id': str})

# norm of x, y, z vector, divide by gravitational acceleration constant, more than 2G -> qc flag
# separate scripts for adding norm/acceleration, flag and filter

data = data[(data['x'].between(*coordinates_span)) & 
            (data['y'].between(*coordinates_span)) &
            (data['z'].between(*coordinates_span))]

data.to_csv(output_file, index=False)
