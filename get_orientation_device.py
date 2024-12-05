import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
paths = pd.read_csv(input_file, dtype=dtype)

summary_list = []
for index, row in paths.iterrows():

    data = pd.read_csv(row.full_path, sep='\t', parse_dates=['timestamp'])
    
    rename = {'value_0': 'x', 'value_1': 'y', 'value_2': 'z'}
    data = data.filter(rename.keys()) \
               .rename(columns=rename)
    
    summary = data[rename.values()].median()
    summary = pd.DataFrame(summary).T
    
    summary['subject_id'] = row.subject_id    
    summary['device_id'] = row.device_id
    summary['test_repeat'] = row.test_repeat
    summary['test_run_type'] = row.test_run_type
    summary['test_run_begin'] = row.test_run_begin
    summary['day_of_study'] = row.day_of_study    
    
    summary_list.append(summary)

summary_tot = pd.concat(summary_list)
summary_tot.to_csv(output_file, index=None, header=True)
