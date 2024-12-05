import pandas as pd

input_duration = snakemake.input.duration[0]
input_acceleration = snakemake.input.acceleration[0]
input_cluster = snakemake.input.cluster[0]
input_steps = snakemake.input.steps[0]
merge_cols = snakemake.params.merge_cols
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
cols = merge_cols + ['qc_pass', 'span_duration_ok']
duration = pd.read_csv(input_duration, dtype=dtype) \
             .filter(cols) \
             .drop_duplicates()

cols = merge_cols + ['acceleration_g_ok']
acceleration = pd.read_csv(input_acceleration, dtype=dtype) \
                 .filter(cols)

cols = merge_cols + ['cluster_ok']
cluster = pd.read_csv(input_cluster, dtype=dtype) \
            .filter(cols)

cols = merge_cols + ['steps_ok']
steps = pd.read_csv(input_steps, dtype=dtype) \
          .filter(cols)

merged = duration.merge(acceleration, on=merge_cols) \
                 .merge(cluster, on=merge_cols) \
                 .merge(steps, on=merge_cols)

merged = merged.rename(columns={'qc_pass': 'phone_not_on_table_ok'})

merged['new_qc_pass'] = merged[['span_duration_ok', 'acceleration_g_ok', 'cluster_ok']] \
                              .sum(axis=1) == 3

merged.to_csv(output_file, index=None, header=True)
