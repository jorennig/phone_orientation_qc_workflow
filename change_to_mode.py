import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-10324-effect-phone-location-balance-test\results\k_means\clusters_combined_path_length.csv'
variable = 'cluster'
variable = 'main_orientation'

data = pd.read_csv(input_file, dtype={'subject_id': str})

cols = ['subject_id']
data['mode_' + variable] = data.groupby(cols)[variable].transform(lambda x: x.value_counts().idxmax())
data['diff_mode_' + variable] = (data['mode_' + variable] != data[variable]).astype(int)
data['diff_total' + variable] = data.groupby(cols)['diff_mode_' + variable].transform('sum')

data.to_csv(output_file, index=None, header=True)
