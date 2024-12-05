import pandas as pd

input_file = snakemake.input.data[0]
grav_acc_threshold = snakemake.params.grav_acc_threshold
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['acceleration_g_ok'] = data['g_force'] <= grav_acc_threshold

data.to_csv(output_file, index=False)
