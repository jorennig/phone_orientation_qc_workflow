import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

input_file = r'\\exports.hps.kau.science.roche.com\pred\rpmda\BN40423-v2\jira\RPMDA-10324-effect-phone-location-balance-test\results\k_means\clusters_combined_changes.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})

data['two_week_period'] = data['day_of_study'] // 14
data = data[data['two_week_period'].between(0, 40)]

group_cols = ['subject_id', 'two_week_period', 'feature_digital']
summary = data.groupby(group_cols) \
              .agg({'change_to_previous': 'sum', 'value_digital': 'median'}) \
              .reset_index() \
              .rename(columns={'change_to_previous': 'orientation_changes'})

summary_orientation = data.groupby(group_cols)['main_orientation'].value_counts() \
                          .reset_index(name='n')
summary_orientation = pd.pivot_table(summary_orientation, index=group_cols,
                                     columns='main_orientation', values='n') \
                        .reset_index() \
                        .fillna(0)

summary_orientation['main_orientation'] = summary_orientation[['x', 'y', 'z']] \
                                          .apply(lambda x: x.argmax(), axis=1)
         
summary_orientation['main_orientation'] = summary_orientation['main_orientation'] \
                                          .replace({0: 'x', 1: 'y', 2: 'z',})

summary_orientation['n_tests'] = summary_orientation[['x', 'y', 'z']] \
                                 .apply(lambda x: x.sum(), axis=1)
summary_orientation['max_orientation'] = summary_orientation[['x', 'y', 'z']] \
                                         .apply(lambda x: x.max(), axis=1)
summary_orientation['percent_main_orientation'] = \
     summary_orientation['max_orientation'] / summary_orientation['n_tests'] * 100

summary_orientation = summary_orientation.filter(['subject_id', 'two_week_period',
                                                  'main_orientation', 'percent_main_orientation'])

summary = summary.merge(summary_orientation, on=group_cols)

summary = pd.melt(summary, id_vars=['subject_id', 'two_week_period', 'feature_digital', 'value_digital'],
                  value_vars=['orientation_changes', 'main_orientation', 'percent_main_orientation'],
                  var_name='feature_orientation', value_name='value_orientation')

summary.to_csv(output_file, index=None, header=True)
