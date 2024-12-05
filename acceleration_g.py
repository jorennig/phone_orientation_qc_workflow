import pandas as pd
import numpy as np

input_file = snakemake.input.data[0]
grav_acc_constant = snakemake.params.grav_acc_constant
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'device_id': str, 'subject_id': str})

data['g_force'] = (np.apply_along_axis(np.linalg.norm, 1, data[['x', 'y', 'z']])) /\
                   grav_acc_constant

data.to_csv(output_file, index=False)
