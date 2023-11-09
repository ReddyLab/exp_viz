use cov_viz_ds::ExperimentFeatureData;

pub fn intersect_coverage_data_features(
    feature_data: Vec<ExperimentFeatureData>,
) -> ExperimentFeatureData {
    if feature_data.len() == 0 {
        return ExperimentFeatureData::default();
    }

    feature_data
        .iter()
        .skip(1)
        .fold(feature_data[0].to_owned(), |acc, feature_data| {
            ExperimentFeatureData {
                sources: acc.sources & &feature_data.sources,
                targets: acc.targets & &feature_data.targets,
            }
        })
}
