import pandas as pd
from glob import glob
import os

configfile: 'config.yaml'

FEATURES = config['features']
DIGITAL_FEATURES = config['digital_features']
CLINICAL_FEATURES = config['clinical_features']
COORDINATES = config['coordinates_dict'].keys()
QC_CHECKS_OLD = config['qc_checks_old']
QC_CHECKS_NEW = config['qc_checks_new']
QC_CHECKS = config['qc_checks']
STEP_FEATURES = config['step_features']

rule format_test_file:
    input: script='scripts/format_test_file.py',
           data=config['balance_file']
    output: temp('data/balance_format.csv')
    params: digital_features=config['digital_features']
    script: 'scripts/format_test_file.py'

rule flag_span_duration:
    input: script='scripts/flag_span_duration.py',
           data=rules.format_test_file.output
    params: span_duration_bounds=config['span_duration_bounds']
    output: 'data/balance_duration_flagged.csv'
    script: 'scripts/flag_span_duration.py'

rule filter_span_duration:
    input: script='scripts/filter_qc.py',
           data=rules.flag_span_duration.output
    params: qc_variable='span_duration_ok'
    output: 'data/balance_duration_filtered.csv'
    script: 'scripts/filter_qc.py'

rule device_data_path:
    input: script='scripts/device_data_path.py',
           data=rules.filter_span_duration.output
    output: 'data/device_data_path.csv'
    params: device_data_path=config['device_data_path']
    script: 'scripts/device_data_path.py'

rule orientation_summary:
    input: script='scripts/get_orientation_device.py',
           data=rules.device_data_path.output
    output: 'data/orientation_summary.csv'
    script: 'scripts/get_orientation_device.py'

rule absolute_values:
    input: script='scripts/absolute_values.py',
           data=rules.orientation_summary.output
    output: 'data/orientation_summary_abs.csv'
    script: 'scripts/absolute_values.py'

rule acceleration_g:
    input: script='scripts/acceleration_g.py',
           data=rules.absolute_values.output
    output: 'data/orientation_summary_g.csv'
    params: grav_acc_constant=config['grav_acc_constant']
    script: 'scripts/acceleration_g.py'

rule flag_acceleration_g:
    input: script='scripts/flag_acceleration_g.py',
           data=rules.acceleration_g.output
    output: 'data/orientation_summary_g_flagged.csv'
    params: grav_acc_threshold=config['grav_acc_threshold']
    script: 'scripts/flag_acceleration_g.py'

rule filter_acceleration_g:
    input: script='scripts/filter_qc.py',
           data=rules.flag_acceleration_g.output
    output: 'data/orientation_summary_g_filtered.csv'
    params: qc_variable='acceleration_g_ok'
    script: 'scripts/filter_qc.py'

rule get_main_orientation:
    input: script='scripts/get_main_orientation.py',
           data=rules.filter_acceleration_g.output
    output: 'data/orientation_summary_main_orientation.csv'
    script: 'scripts/get_main_orientation.py'

rule change_to_previous_main_orientation:
    input: script='scripts/change_to_previous.py',
           data=rules.get_main_orientation.output,
    output: 'data/orientation_summary_change_to_previous.csv'
    params: variable='main_orientation'
    script: 'scripts/change_to_previous.py'

rule histogram:
    input: script='scripts/histogram.py',
           data=rules.change_to_previous_main_orientation.output
    output: 'results/histogram_{dimension}_coordinate.png'
    script: 'scripts/histogram.py'

rule pca_coordinates:
    input: script='scripts/pca_coordinates.py',
           data=rules.change_to_previous_main_orientation.output
    output: data='data/orientation_summary_pca.csv',
            variance_explained='results/pca_variance_explained.csv'
    script: 'scripts/pca_coordinates.py'

rule bar_plot_variance_explained:
    input: script='scripts/bar_plot.py',
           data=rules.pca_coordinates.output.variance_explained
    output: 'results/bar_plot_variance_explained.png'
    params: x='pc', y='variance_explained'
    script: 'scripts/bar_plot.py'

rule k_means_coordinates:
    input: script='scripts/k_means_coordinates.py',
           data=rules.pca_coordinates.output.data
    output: 'results/k_means/{n_clust}_clusters.csv'
    params: kmeans_args=config['kmeans_args']
    script: 'scripts/k_means_coordinates.py'

rule elbow_plot:
    input: script='scripts/elbow_plot.py',
           data=rules.pca_coordinates.output.data
    output: 'results/k_means/elbow_plot.png'
    params: kmeans_args=config['kmeans_args'],
            max_clust=config['max_clust']
    script: 'scripts/elbow_plot.py'

rule scatterplot_pca:
    input: script='scripts/scatterplot_pca.py',
           data=rules.k_means_coordinates.output
    output: 'results/k_means/scatterplot_pca_{n_clust}_clusters.png'
    script: 'scripts/scatterplot_pca.py'

rule density_2d_clusters:
    input: script='scripts/density_2d.py',
           data=rules.k_means_coordinates.output
    output: 'results/k_means/density_2d_{n_clust}_clusters_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='cluster'
    script: 'scripts/density_2d.py'

rule bin_2d_clusters:
    input: script='scripts/bin_2d.py',
           data=rules.k_means_coordinates.output
    output: 'results/k_means/bin_2d_{n_clust}_clusters_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='cluster'
    script: 'scripts/bin_2d.py'

rule density_2d_main_orientation:
    input: script='scripts/density_2d.py',
           data=rules.get_main_orientation.output
    output: 'results/main_orientation/density_2d_main_orientation_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='main_orientation'
    script: 'scripts/density_2d.py'

rule bin_2d_main_orientation:
    input: script='scripts/bin_2d.py',
           data=rules.get_main_orientation.output
    output: 'results/main_orientation/bin_2d_main_orientation_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='main_orientation'
    script: 'scripts/bin_2d.py'

rule combine_clusters:
    input: script='scripts/combine_clusters.py',
           data='results/k_means/4_clusters.csv'
    output: 'results/k_means/clusters_combined.csv'
    params: combine_clusters=config['combine_clusters']
    script: 'scripts/combine_clusters.py'

rule flag_cluster:
    input: script='scripts/flag_cluster.py',
           data=rules.combine_clusters.output
    params: cluster_flag=config['cluster_flag']
    output: 'results/k_means/clusters_combined_flagged.csv'
    script: 'scripts/flag_cluster.py'

rule get_main_cluster:
    input: script='scripts/get_main_cluster.py',
           data=rules.flag_cluster.output
    params: cluster_flag=config['cluster_flag']
    output: 'results/k_means/clusters_combined_main_cluster.csv'
    script: 'scripts/get_main_cluster.py'

rule change_to_previous_cluster:
    input: script='scripts/change_to_previous.py',
           data=rules.get_main_cluster.output,
    output: 'results/k_means/clusters_combined_change_to_previous.csv'
    params: variable='cluster'
    script: 'scripts/change_to_previous.py'

rule scatterplot_mismatch_confusion_matrix:
    input: script='scripts/scatterplot_mismatch_confusion_matrix.py',
           data=rules.change_to_previous_cluster.output
    output: 'results/k_means/scatterplot_mismatch_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            cluster_mismatch=config['cluster_mismatch'],
            orientation_mismatch=config['orientation_mismatch'],
    script: 'scripts/scatterplot_mismatch_confusion_matrix.py'

rule scatterplot_pca_cluster_combined:
    input: script='scripts/scatterplot_pca.py',
           data=rules.change_to_previous_cluster.output
    output: 'results/k_means/scatterplot_pca_clusters_combined.png'
    script: 'scripts/scatterplot_pca.py'

rule density_2d_clusters_combined:
    input: script='scripts/density_2d.py',
           data=rules.change_to_previous_cluster.output
    output: 'results/k_means/density_2d_clusters_combined_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='cluster'
    script: 'scripts/density_2d.py'

rule bin_2d_clusters_combined:
    input: script='scripts/bin_2d.py',
           data=rules.change_to_previous_cluster.output
    output: 'results/k_means/bin_2d_clusters_combined_{coords}.png'
    params: coordinates_dict=config['coordinates_dict'],
            color='cluster'
    script: 'scripts/bin_2d.py'

rule summary_study:
    input: script='scripts/summary_study.py',
           data=rules.change_to_previous_cluster.output
    output: table='results/summary_study.csv'
    script: 'scripts/summary_study.py'

rule bar_plot_summary_study:
    input: script='scripts/bar_plot.py',
           data=rules.summary_study.output.table
    output: 'results/bar_plot_summary_study.png'
    params: x='cluster', y='n'
    script: 'scripts/bar_plot.py'

rule summary_subject:
    input: script='scripts/summary_subject.py',
           data=rules.change_to_previous_cluster.output
    output: 'results/summary_subject.csv'
    params: features=config['features']
    script: 'scripts/summary_subject.py'

rule summary_subject_n_tests_threshold_flag:
    input: script='scripts/n_tests_threshold_flag.py',
           data=rules.summary_subject.output
    output: 'results/summary_subject_n_tests_flagged.csv'
    params: quantile=config['quantile_threshold']
    script: 'scripts/n_tests_threshold_flag.py'

rule summary_subject_n_tests_threshold_filter:
    input: script='scripts/filter_qc.py',
           data=rules.summary_subject_n_tests_threshold_flag.output
    output: 'results/summary_subject_n_tests_filtered.csv'
    params: qc_variable='n_tests_ok'
    script: 'scripts/filter_qc.py'

rule merge_clusters_features:
    input: script='scripts/merge_features.py',
           data=rules.change_to_previous_cluster.output,
           test=rules.filter_span_duration.output
    output: 'results/k_means/clusters_combined_digital_features.csv'
    params: merge_cols=config['merge_cols']
    script: 'scripts/merge_features.py'

rule summary_cluster_digital:
    input: script='scripts/summary_cluster.py',
           data=rules.merge_clusters_features.output,
    output: 'results/k_means/analysis_digital/summary_cluster_digital.csv'
    params: feature='feature_digital',
            value='value_digital',
            cluster_col='cluster'
    script: 'scripts/summary_cluster.py'

rule histogram_cluster_digital:
    input: script='scripts/histogram_cluster.py',
           data=rules.merge_clusters_features.output,
    output: 'results/k_means/analysis_digital/histogram_cluster_{feature_test}.png'
    params: feature='feature_digital',
            value='value_digital',
            cluster_col='cluster'
    script: 'scripts/histogram_cluster.py'

rule box_plot_cluster_digital:
    input: script='scripts/box_plot_cluster.py',
           data=rules.merge_clusters_features.output,
    output: 'results/k_means/analysis_digital/box_plot_cluster_digital.png'
    params: feature='feature_digital',
            value='value_digital',
            cluster_col='cluster',
            alpha=0.05
    script: 'scripts/box_plot_cluster.py'

rule kruskal_wallis_cluster_digital:
    input: script='scripts/kruskal_wallis_cluster.py',
           data=rules.merge_clusters_features.output,
    output: 'results/k_means/analysis_digital/kruskal_wallis_cluster_digital.csv'
    params: feature='feature_digital',
            value='value_digital',
            cluster_col='cluster'
    script: 'scripts/kruskal_wallis_cluster.py'

rule u_test_cluster_digital:
    input: script='scripts/u_test_cluster.py',
           data=rules.merge_clusters_features.output,
    output: temp('results/k_means/analysis_digital/u_test_cluster_digital.csv')
    params: feature='feature_digital',
            value='value_digital',
            cluster_col='cluster'
    script: 'scripts/u_test_cluster.py'

rule effect_size_cluster_digital:
    input: script='scripts/effect_size.py',
           data=rules.merge_clusters_features.output,
    output: temp('results/k_means/analysis_digital/effect_size_cluster_digital.csv')
    params: feature='feature_digital',
            value='value_digital',
            cluster='cluster'
    script: 'scripts/effect_size.py'

rule merge_u_test_effect_size_digital:
    input: script='scripts/merge_u_test_effect_size.py',
           u_test=rules.u_test_cluster_digital.output,
           effect_size=rules.effect_size_cluster_digital.output,
    output: 'results/k_means/analysis_digital/u_test_effect_size_cluster_digital.csv'
    params: feature='feature_digital',
    script: 'scripts/merge_u_test_effect_size.py'

rule confusion_matrix_clusters:
    input: script='scripts/confusion_matrix_clusters.py',
           data=rules.change_to_previous_cluster.output,
    output: 'results/k_means/confusion_matrix_clusters.csv'
    script: 'scripts/confusion_matrix_clusters.py'

rule plot_confusion_matrix_clusters:
    input: script='scripts/plot_confusion_matrix_clusters.py',
           data=rules.confusion_matrix_clusters.output,
    output: 'results/k_means/confusion_matrix_clusters.png'
    script: 'scripts/plot_confusion_matrix_clusters.py'

rule lmem_orientation:
    input: script='scripts/lmem_orientation.R',
           data=rules.merge_clusters_features.output,
    output: 'results/lmem_orientation.csv'
    script: 'scripts/lmem_orientation.R'

rule select_clinical:
    input: script='scripts/select_clinical.py',
           data=config['clinical']
    output: 'data/clinical_data.csv'
    params: clinical_features=config['clinical_features'],
    script: 'scripts/select_clinical.py'

rule merge_clinical:
    input: script='scripts/merge_clinical.py',
           data=rules.summary_subject_n_tests_threshold_filter.output,
           clinical=rules.select_clinical.output,
    output: 'results/summary_subject_clinical.csv'
    script: 'scripts/merge_clinical.py'

rule correlation_clinical:
    input: script='scripts/correlation.py',
           data=rules.merge_clinical.output,
    output: 'results/correlation_clinical_balance.csv'
    params: group_cols=['feature_clinical', 'feature_orientation', 'qc_type'],
            value_cols=['value_clinical', 'value_orientation']
    script: 'scripts/correlation.py'
 
rule cluster_clinical:
    input: script='scripts/cluster_clinical.py',
           data=rules.merge_clinical.output,
    output: 'results/k_means/cluster_clinical.csv'
    script: 'scripts/cluster_clinical.py'

rule summary_cluster_clinical:
    input: script='scripts/summary_cluster.py',
           data=rules.cluster_clinical.output,
    output: 'results/k_means/analysis_clinical/summary_cluster_clinical.csv'
    params: feature='feature_clinical',
            value='value_clinical',
            cluster_col='main_cluster'
    script: 'scripts/summary_cluster.py'

rule histogram_cluster_clinical:
    input: script='scripts/histogram_cluster.py',
           data=rules.cluster_clinical.output,
    output: 'results/k_means/analysis_clinical/histogram_cluster_{feature_test}.png'
    params: feature='feature_clinical',
            value='value_clinical', 
            cluster_col='main_cluster'            
    script: 'scripts/histogram_cluster.py'

rule box_plot_cluster_clinical:
    input: script='scripts/box_plot_cluster.py',
           data=rules.cluster_clinical.output,
    output: 'results/k_means/analysis_clinical/box_plot_cluster_clinical.png'
    params: feature='feature_clinical',
            value='value_clinical',
            cluster_col='main_cluster',
            alpha=0.4
    script: 'scripts/box_plot_cluster.py'

rule kruskal_wallis_cluster_clinical:
    input: script='scripts/kruskal_wallis_cluster.py',
           data=rules.cluster_clinical.output,
    output: 'results/k_means/analysis_clinical/kruskal_wallis_cluster_clinical.csv'
    params: feature='feature_clinical',
            value='value_clinical',
            cluster_col='main_cluster'
    script: 'scripts/kruskal_wallis_cluster.py'

rule u_test_cluster_clinical:
    input: script='scripts/u_test_cluster.py',
           data=rules.cluster_clinical.output,
    output: temp('results/k_means/analysis_clinical/u_test_cluster_clinical.csv')
    params: feature='feature_clinical',
            value='value_clinical',
            cluster_col='main_cluster'
    script: 'scripts/u_test_cluster.py'

rule effect_size_cluster_clinical:
    input: script='scripts/effect_size.py',
           data=rules.cluster_clinical.output,
    output: temp('results/k_means/analysis_clinical/effect_size_cluster_clinical.csv')
    params: feature='feature_clinical',
            value='value_clinical',
            cluster='main_cluster'
    script: 'scripts/effect_size.py'

rule merge_u_test_effect_size_clinical:
    input: script='scripts/merge_u_test_effect_size.py',
           u_test=rules.u_test_cluster_clinical.output,
           effect_size=rules.effect_size_cluster_clinical.output,
    output: 'results/k_means/analysis_clinical/u_test_effect_size_cluster_clinical.csv'
    params: feature='feature_clinical',
    script: 'scripts/merge_u_test_effect_size.py'

rule scatterplot_clinical:
    input: script='scripts/scatterplot_clinical.py',
           data=rules.merge_clinical.output,
           results=rules.correlation_clinical.output
    output: 'results/scatterplots_clinical/{feature}_clinical.png'
    params: feature_type='feature_orientation',
            value_type='value_orientation',
            p_threshold=config['p_threshold'],
            merge_cols=['feature_clinical', 'feature_orientation', 'qc_type'],
    script: 'scripts/scatterplot_clinical.py'

rule format_steps_qc:
    input: script='scripts/format_steps_qc.py',
           data=config['balance_file']
    output: 'data/balance_steps_qc.csv'
    params: steps_qc=config['steps_qc']
    script: 'scripts/format_steps_qc.py'

rule merge_qc:
    input: script='scripts/merge_qc.py',
           duration=rules.flag_span_duration.output,
           acceleration=rules.flag_acceleration_g.output,
           cluster=rules.change_to_previous_cluster.output,
           steps=rules.format_steps_qc.output,
    output: 'results/merge_qc.csv'
    params: merge_cols=config['merge_cols']
    script: 'scripts/merge_qc.py'

rule qc_pass_rate:
    input: script='scripts/qc_pass_rate.py',
           data=rules.merge_qc.output
    output: 'data/qc_pass_rate.csv'
    script: 'scripts/qc_pass_rate.py'

rule histogram_n_tests:
    input: script='scripts/histogram_variable.py',
           data=rules.qc_pass_rate.output
    output: 'results/histogram_n_tests.png'
    params: variable='n_tests'
    script: 'scripts/histogram_variable.py'

rule n_tests_threshold_flag:
    input: script='scripts/n_tests_threshold_flag.py',
           data=rules.qc_pass_rate.output
    output: 'data/qc_pass_rate_n_tests_flagged.csv'
    params: quantile=config['quantile_threshold']
    script: 'scripts/n_tests_threshold_flag.py'

rule n_tests_threshold_filter:
    input: script='scripts/filter_qc.py',
           data=rules.n_tests_threshold_flag.output
    output: 'data/qc_pass_rate_n_tests_filtered.csv'
    params: qc_variable='n_tests_ok'
    script: 'scripts/filter_qc.py'

rule merge_qc_pass_rate_clinical:
    input: script='scripts/merge_clinical.py',
           data=rules.n_tests_threshold_filter.output,
           clinical=rules.select_clinical.output,
    output: 'data/qc_pass_rate_clinical.csv'
    script: 'scripts/merge_clinical.py'

rule correlation_qc_pass_rate_clinical:
    input: script='scripts/correlation.py',
           data=rules.merge_qc_pass_rate_clinical.output,
    output: 'results/correlation_qc_pass_rate_clinical.csv'
    params: group_cols=['feature_clinical', 'feature_qc', 'qc_type'],
            value_cols=['value_clinical', 'qc_pass_rate_percent']
    script: 'scripts/correlation.py'

rule scatterplot_qc_pass_rate_clinical:
    input: script='scripts/scatterplot_clinical.py',
           data=rules.merge_qc_pass_rate_clinical.output,
           results=rules.correlation_qc_pass_rate_clinical.output
    output: 'results/scatterplots_qc_pass_rate_clinical/{feature}_clinical.png'
    params: feature_type='feature_qc',
            value_type='qc_pass_rate_percent',
            p_threshold=config['p_threshold'],
            merge_cols=['feature_clinical', 'feature_qc', 'qc_type'],
    script: 'scripts/scatterplot_clinical.py'

rule confusion_matrix_qc:
    input: script='scripts/confusion_matrix_qc.py',
           data=rules.merge_qc.output,
    output: 'results/confusion_matrix_qc/confusion_matrix-{qc_check_old}-{qc_check_new}.csv'
    params: merge_cols=config['merge_cols']
    script: 'scripts/confusion_matrix_qc.py'

rule plot_confusion_matrix_qc:
    input: script='scripts/plot_confusion_matrix_qc.py',
           data=rules.confusion_matrix_qc.output,
    output: 'results/confusion_matrix_qc/confusion_matrix-{qc_check_old}-{qc_check_new}.png'
    script: 'scripts/plot_confusion_matrix_qc.py'

rule merge_qc_digital:
    input: script='scripts/merge_qc_digital.py',
           qc=rules.merge_qc.output,
           digital=rules.flag_span_duration.output,
    output: 'results/merge_qc_digital.csv'
    params: merge_cols=config['merge_cols'],
    script: 'scripts/merge_qc_digital.py'

rule u_test_qc_digital:
    input: script='scripts/u_test_qc_digital.py',
           data=rules.merge_qc_digital.output,
    output: 'results/u_test_qc_digital.csv'
    script: 'scripts/u_test_qc_digital.py'

rule box_plot_qc_digital:
    input: script='scripts/box_plot_qc_digital.py',
           data=rules.merge_qc_digital.output,
    output: 'results/box_plot_qc_digital.png'
    script: 'scripts/box_plot_qc_digital.py'

rule filter_qc_digital:
    input: script='scripts/filter_qc.py',
           data=rules.merge_qc_digital.output,
    output: 'results/merge_qc_digital_filtered.csv'
    params: qc_variable='qc_value'
    script: 'scripts/filter_qc.py'

rule aggregate_qc_digital:
    input: script='scripts/aggregate_qc_digital.py',
           data=rules.filter_qc_digital.output,
    output: 'results/merge_qc_digital_aggregated.csv'
    script: 'scripts/aggregate_qc_digital.py'

rule merge_clinical_qc_digital:
    input: script='scripts/merge_clinical.py',
           data=rules.aggregate_qc_digital.output,
           clinical=rules.select_clinical.output,
    output: 'results/merge_qc_digital_clinical.csv'
    script: 'scripts/merge_clinical.py'

rule correlation_clinical_qc_digital:
    input: script='scripts/correlation.py',
           data=rules.merge_clinical_qc_digital.output,
    output: 'results/correlation_clinical_qc_digital.csv'
    params: group_cols=['feature_clinical', 'feature_digital', 'qc_type'],
            value_cols=['value_clinical', 'value_digital']
    script: 'scripts/correlation.py'

rule scatterplot_clinical_qc_digital:
    input: script='scripts/scatterplot_clinical.py',
           data=rules.merge_clinical_qc_digital.output,
           results=rules.correlation_clinical_qc_digital.output
    output: 'results/scatterplots_clinical_qc_digital/{feature}_clinical.png'
    params: feature_type='feature_digital',
            value_type='value_digital',
            p_threshold=config['p_threshold'],
            merge_cols=['feature_clinical', 'feature_digital', 'qc_type'],
    script: 'scripts/scatterplot_clinical.py'

rule compare_qc:
    input: script='scripts/compare_qc.py',
           data=rules.aggregate_qc_digital.output,
    output: 'results/compare_qc.csv'
    script: 'scripts/compare_qc.py'

rule histogram_compare_qc:
    input: script='scripts/histogram_variable.py',
           data=rules.compare_qc.output
    output: 'results/histogram_compare_qc.png'
    params: variable='diff_value'
    script: 'scripts/histogram_variable.py'

rule prepare_steps:
    input: script='scripts/prepare_steps.py',
           data=rules.format_steps_qc.output,
    output: 'data/balance_steps_aggregated.csv'
    script: 'scripts/prepare_steps.py'

rule merge_steps_clinical:
    input: script='scripts/merge_clinical.py',
           data=rules.prepare_steps.output,
           clinical=rules.select_clinical.output,
    output: 'data/balance_steps_clinical.csv'
    script: 'scripts/merge_clinical.py'

rule correlation_steps_clinical:
    input: script='scripts/correlation.py',
           data=rules.merge_steps_clinical.output,
    output: 'results/correlation_steps_clinical.csv'
    params: group_cols=['feature_clinical', 'feature_digital', 'qc_type'],
            value_cols=['value_clinical', 'value_digital']
    script: 'scripts/correlation.py'

rule scatterplot_steps_digital:
    input: script='scripts/scatterplot_clinical.py',
           data=rules.merge_steps_clinical.output,
           results=rules.correlation_steps_clinical.output
    output: 'results/scatterplots_clinical_steps/{feature}_clinical.png'
    params: feature_type='feature_digital',
            value_type='value_digital',
            p_threshold=config['p_threshold'],
            merge_cols=['feature_clinical', 'feature_digital', 'qc_type'],
    script: 'scripts/scatterplot_clinical.py'

rule prepare_cluster_data:
    input: script='scripts/prepare_cluster_data.py',
           data=rules.summary_subject_n_tests_threshold_filter.output,
    output: 'results/summary_subject_main_cluster.csv'
    script: 'scripts/prepare_cluster_data.py'

rule kruskal_wallis_percent_main_cluster:
    input: script='scripts/kruskal_wallis_percent_main_cluster.py',
           data=rules.prepare_cluster_data.output,
    output: 'results/k_means/kruskal_wallis_percent_main_cluster.csv'
    script: 'scripts/kruskal_wallis_percent_main_cluster.py'

rule box_plot_percent_main_cluster:
    input: script='scripts/box_plot_percent_main_cluster.py',
           data=rules.prepare_cluster_data.output,
    output: 'results/k_means/box_plot_percent_main_cluster.png'
    script: 'scripts/box_plot_percent_main_cluster.py'

rule u_test_percent_main_cluster:
    input: script='scripts/u_test_percent_main_cluster.py',
           data=rules.prepare_cluster_data.output,
    output: 'results/k_means/u_test_percent_main_cluster.csv'
    script: 'scripts/u_test_percent_main_cluster.py'


rule all:
    input: expand('results/histogram_{dimension}_coordinate.png', dimension=config['dimensions']),
           rules.bar_plot_variance_explained.output,
           expand('results/k_means/{n_clust}_clusters.csv', n_clust=config['n_clusters']),
           rules.elbow_plot.output,
           expand('results/k_means/scatterplot_pca_{n_clust}_clusters.png', n_clust=config['n_clusters']),
           expand('results/k_means/density_2d_{n_clust}_clusters_{coords}.png', n_clust=config['n_clusters'], coords=COORDINATES),
           expand('results/k_means/bin_2d_{n_clust}_clusters_{coords}.png', n_clust=config['n_clusters'], coords=COORDINATES),
           expand('results/main_orientation/density_2d_main_orientation_{coords}.png', coords=COORDINATES),
           expand('results/main_orientation/bin_2d_main_orientation_{coords}.png', coords=COORDINATES),
           rules.change_to_previous_cluster.output,
           expand('results/k_means/scatterplot_mismatch_{coords}.png', coords=COORDINATES),
           rules.scatterplot_pca_cluster_combined.output,
           expand('results/k_means/density_2d_clusters_combined_{coords}.png', coords=COORDINATES),
           expand('results/k_means/bin_2d_clusters_combined_{coords}.png', coords=COORDINATES),
           rules.bar_plot_summary_study.output,
           rules.summary_cluster_digital.output,
           expand('results/k_means/analysis_digital/histogram_cluster_{feature_test}.png', feature_test=DIGITAL_FEATURES),
           rules.box_plot_cluster_digital.output,
           rules.kruskal_wallis_cluster_digital.output,
           rules.merge_u_test_effect_size_digital.output,
           rules.plot_confusion_matrix_clusters.output,
           rules.lmem_orientation.output,
           rules.summary_cluster_clinical.output,
           expand('results/k_means/analysis_clinical/histogram_cluster_{feature_test}.png', feature_test=CLINICAL_FEATURES),
           rules.box_plot_cluster_clinical.output,
           rules.kruskal_wallis_cluster_clinical.output,
           rules.merge_u_test_effect_size_clinical.output,
           expand('results/scatterplots_clinical/{feature}_clinical.png', feature=FEATURES),
           rules.format_steps_qc.output,
           rules.merge_qc.output,
           rules.histogram_n_tests.output,           
           rules.correlation_qc_pass_rate_clinical.output,
           expand('results/scatterplots_qc_pass_rate_clinical/{feature}_clinical.png', feature=QC_CHECKS),
           expand('results/confusion_matrix_qc/confusion_matrix-{qc_check_old}-{qc_check_new}.png', qc_check_old=QC_CHECKS_OLD, qc_check_new=QC_CHECKS_NEW),
           rules.u_test_qc_digital.output,
           rules.box_plot_qc_digital.output,
           rules.correlation_clinical_qc_digital.output,
           expand('results/scatterplots_clinical_qc_digital/{feature}_clinical.png', feature=DIGITAL_FEATURES),
           rules.histogram_compare_qc.output,
           rules.correlation_steps_clinical.output,
           expand('results/scatterplots_clinical_steps/{feature}_clinical.png', feature=STEP_FEATURES),
           rules.kruskal_wallis_percent_main_cluster.output,
           rules.box_plot_percent_main_cluster.output,
           rules.u_test_percent_main_cluster.output,
           