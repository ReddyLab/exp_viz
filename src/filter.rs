use rayon::prelude::*;
use roaring::RoaringTreemap;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::filter_data_structures::*;
use cov_viz_ds::{
    BucketLoc, CoverageData, DbID, ExperimentFeatureData, FacetRange, FacetRange64, ObservationData,
};

#[derive(Debug, Clone)]
struct BucketData {
    feature_ids: RoaringTreemap,
    associated_features: RoaringTreemap,
    min_effect: f32,
    max_effect: f32,
    min_sig: f64,
    max_sig: f64,
}

fn is_disjoint(a: &Vec<DbID>, b: &Vec<DbID>) -> bool {
    for val_a in a {
        for val_b in b {
            if val_a == val_b {
                return false;
            }
        }
    }

    true
}

fn add_data_to_bucket(
    id: DbID,
    associated_feature: Option<DbID>,
    obs_sig: f64,
    effect_size: f32,
    buckets: &mut FxHashMap<BucketLoc, BucketData>,
    bucket_locs: &FxHashMap<DbID, BucketLoc>,
) {
    let feature_loc = bucket_locs.get(&id);
    if feature_loc.is_none() {
        return;
    }

    buckets
        .entry(*feature_loc.unwrap())
        .and_modify(|bucket_data| {
            bucket_data.feature_ids.insert(id);
            if let Some(af) = associated_feature {
                bucket_data.associated_features.insert(af);
            }
            bucket_data.min_effect = effect_size.min(bucket_data.min_effect);
            bucket_data.max_effect = effect_size.max(bucket_data.max_effect);

            bucket_data.min_sig = obs_sig.min(bucket_data.min_sig);
            bucket_data.max_sig = obs_sig.max(bucket_data.max_sig);
        })
        .or_insert(BucketData {
            feature_ids: RoaringTreemap::from([id]),
            associated_features: if let Some(af) = associated_feature {
                RoaringTreemap::from([af])
            } else {
                RoaringTreemap::new()
            },
            min_effect: effect_size,
            max_effect: effect_size,
            min_sig: obs_sig,
            max_sig: obs_sig,
        });
}

fn update_buckets(
    observation: &ObservationData,
    source_buckets: &mut FxHashMap<BucketLoc, BucketData>,
    target_buckets: &mut FxHashMap<BucketLoc, BucketData>,
    features: &FxHashMap<DbID, BucketLoc>,
) {
    add_data_to_bucket(
        observation.source_id,
        observation.target_id,
        observation.neg_log_significance,
        observation.effect_size,
        source_buckets,
        features,
    );

    if let Some(id) = observation.target_id {
        add_data_to_bucket(
            id,
            Some(observation.source_id),
            observation.neg_log_significance,
            observation.effect_size,
            target_buckets,
            features,
        );
    };
}

fn update_bucket_map(
    bucket_list1: &mut FxHashMap<BucketLoc, BucketData>,
    bucket_list2: &FxHashMap<BucketLoc, BucketData>,
) {
    for (loc, data2) in bucket_list2 {
        bucket_list1
            .entry(*loc)
            .and_modify(|bucket_data| {
                bucket_data.feature_ids.extend(&data2.feature_ids);
                bucket_data
                    .associated_features
                    .extend(&data2.associated_features);
                bucket_data.min_effect = data2.min_effect.min(bucket_data.min_effect);
                bucket_data.max_effect = data2.max_effect.max(bucket_data.max_effect);

                bucket_data.min_sig = data2.min_sig.min(bucket_data.min_sig);
                bucket_data.max_sig = data2.max_sig.max(bucket_data.max_sig);
            })
            .or_insert(data2.clone());
    }
}

fn gen_filtered_data(
    buckets: FxHashMap<BucketLoc, BucketData>,
    chrom: Option<u8>,
    feature_count: &mut RoaringTreemap,
    min_effect: &mut f32,
    max_effect: &mut f32,
    min_sig: &mut f64,
    max_sig: &mut f64,
    intervals: &mut Vec<&mut Vec<FilteredBucket>>,
    bucket_size: u32,
    features: &FxHashMap<DbID, BucketLoc>,
) {
    let mut ordered_buckets: Vec<_> = buckets
        .into_iter()
        .filter(|(bucket_loc, _)| chrom.is_none() || bucket_loc.chrom == chrom.unwrap())
        .collect();
    ordered_buckets.sort_by(|(loc1, _), (loc2, _)| loc1.cmp(loc2));
    for (bucket_loc, bucket_data) in ordered_buckets {
        feature_count.extend(&bucket_data.feature_ids);
        *min_effect = min_effect.min(bucket_data.min_effect);
        *max_effect = max_effect.max(bucket_data.max_effect);
        *min_sig = min_sig.min(bucket_data.min_sig);
        *max_sig = max_sig.max(bucket_data.max_sig);

        let chrom = if chrom.is_none() {
            bucket_loc.chrom as usize
        } else {
            0
        };

        intervals[chrom].push(FilteredBucket {
            start: bucket_size * bucket_loc.idx + 1,
            count: bucket_data.feature_ids.len() as usize,
            // buckets are stored as a list where the chromosome indexes and bucket indexes alternate.
            // This cuts down on how much data get sent over the wire.
            associated_buckets: bucket_data
                .associated_features
                .iter()
                .fold(FxHashSet::<&BucketLoc>::default(), |mut acc, id| {
                    if let Some(bucket) = features.get(&id) {
                        acc.insert(bucket);
                    }
                    acc
                })
                .iter()
                .fold(Vec::new(), |mut acc, bucket| {
                    acc.push(bucket.chrom as u32);
                    acc.push(bucket.idx);
                    acc
                }),
            max_log10_sig: bucket_data.max_sig,
            max_abs_effect: if bucket_data.max_effect > bucket_data.min_effect.abs() {
                bucket_data.max_effect
            } else {
                bucket_data.min_effect
            },
        })
    }
}

pub fn filter_coverage_data(
    filters: &Filter,
    data: &CoverageData,
    included_features: Option<&ExperimentFeatureData>,
) -> FilteredData {
    let bucket_size = data.bucket_size;
    let feature_buckets = &data.feature_buckets;

    //
    // Get Numeric Facet Info
    //
    let skip_cont_facet_check = filters.numeric_intervals.is_none();

    let effect_size_interval = match &filters.numeric_intervals {
        Some(c) => FacetRange(c.effect.0, c.effect.1),
        None => data
            .facets
            .iter()
            .find(|f| f.name == "Effect Size")
            .unwrap()
            .range
            .unwrap(),
    };
    let sig_interval = match &filters.numeric_intervals {
        Some(c) => FacetRange64(c.sig.0, c.sig.1),
        None => data
            .facets
            .iter()
            .find(|f| f.name == "Significance")
            .unwrap()
            .range64
            .unwrap(),
    };

    //
    // Get Categorical Facet Info
    //

    // All categorical facet value database ids that are used in this data set.
    // This may not be all possible facet values.
    let mut all_coverage_data_cat_facets: FxHashSet<DbID> = FxHashSet::default();
    for facet in data.facets.iter() {
        if let Some(facet_values) = &facet.values {
            facet_values.keys().for_each(|key| {
                all_coverage_data_cat_facets.insert(*key);
            });
        }
    }

    // Categorical facet value database ids for that are filtered on, not including
    // facet values that aren't used in this data set.
    let coverage_data_cat_facets: FxHashSet<DbID> = all_coverage_data_cat_facets
        .intersection(&filters.categorical_facets)
        .cloned()
        .collect();

    let skip_cat_facet_check = coverage_data_cat_facets.is_empty();

    // All categorical facet value database ids that are used in this data set,
    // divided up by facet
    let mut facet_ids: Vec<FxHashSet<DbID>> = Vec::new();
    for facet in data
        .facets
        .iter()
        .filter(|f| f.facet_type == "FacetType.CATEGORICAL")
    {
        facet_ids.push(FxHashSet::from_iter(
            facet.values.as_ref().unwrap().keys().cloned(),
        ));
    }

    // Facet id sets that have values being filtered on
    let f_with_selections: Vec<FxHashSet<DbID>> = facet_ids
        .into_iter()
        .filter(|f| !f.is_disjoint(&coverage_data_cat_facets))
        .collect();

    // Facet id sets that have values being filtered on and only have the filtered values included in the set
    let selected_f: Vec<Vec<DbID>> = f_with_selections
        .iter()
        .map(|f| (f & &coverage_data_cat_facets).iter().cloned().collect())
        .collect();

    // println!("{:?}", filters.categorical_facets); // all filtered facet values
    // println!("{:?}", all_coverage_data_cat_facets); // all facet values used in data
    // println!("{:?}", coverage_data_cat_facets); // interesection of the above two
    // println!("{:?}", f_with_selections); // Facet id sets that have values being filtered on
    // println!("{:?}", selected_f); // the above, but only including values in the filter

    let direction_facet = data.facets.iter().find(|f| f.name == "Direction").unwrap();
    let direction_facet_values: FxHashSet<DbID> = direction_facet
        .values
        .as_ref()
        .unwrap()
        .iter()
        .map(|fv| *fv.0)
        .collect();
    let nonsignificant_facet_value = direction_facet
        .values
        .as_ref()
        .unwrap()
        .iter()
        .find(|(_, fv_name)| *fv_name == "Non-significant");

    // Skip filtering (i.e., drop completely) non-significant observations IF
    // * at least one direction facet value is checked
    // * and the non-significant facet value isn't checked
    // * There are no non-significant observations
    let skip_nonsignificants = if let Some(nfv) = nonsignificant_facet_value {
        f_with_selections.contains(&direction_facet_values)
            && !coverage_data_cat_facets.contains(&nfv.0)
    } else {
        true
    };

    //
    // Filter Observations
    //

    let empty_vec = Vec::<ObservationData>::new();
    let observations = if skip_nonsignificants {
        data.significant_observations
            .par_iter()
            .chain(empty_vec.par_iter())
    } else {
        data.significant_observations
            .par_iter()
            .chain(data.nonsignificant_observations.par_iter())
    };

    let filtered_observations: Vec<&ObservationData> =
        if skip_cont_facet_check && f_with_selections.is_empty() {
            if let Some(included_features) = included_features {
                observations
                    .filter(|o| -> bool {
                        if let Some(target_id) = o.target_id {
                            included_features.targets.contains(target_id)
                        } else {
                            false
                        }
                    })
                    .collect()
            } else {
                observations.collect()
            }
        } else {
            let filtered_observations = observations.filter(|observation| -> bool {
                if skip_cat_facet_check
                    || selected_f
                        .iter()
                        .all(|f| !is_disjoint(&observation.facet_value_ids, f))
                {
                    if skip_cont_facet_check
                        || (observation.effect_size >= effect_size_interval.0
                            && observation.effect_size <= effect_size_interval.1
                            && observation.neg_log_significance >= sig_interval.0
                            && observation.neg_log_significance <= sig_interval.1)
                    {
                        return true;
                    }
                }

                false
            });

            if let Some(included_features) = included_features {
                filtered_observations
                    .filter(|o| -> bool {
                        if let Some(target_id) = o.target_id {
                            included_features.targets.contains(target_id)
                        } else {
                            false
                        }
                    })
                    .collect()
            } else {
                filtered_observations.collect()
            }
        };

    //
    // Build intermediate bucket data
    //

    // How many "chunks" of work we want to split building buckets into, since it can't be fully parallelized
    // due to using a shared data structure
    let p_count = if let Ok(p_count) = std::thread::available_parallelism() {
        p_count.get()
    } else {
        6
    };

    // Merge filtered observations into an intermediate set of data structures
    // that will then be turned into FilteredData
    let observation_chunks =
        filtered_observations.par_chunks(1.max(filtered_observations.len() / p_count));
    let filter_results: Vec<(
        RoaringTreemap,
        FxHashMap<BucketLoc, BucketData>,
        FxHashMap<BucketLoc, BucketData>,
    )> = observation_chunks
        .map(|chunk| {
            let mut reos = RoaringTreemap::new();
            let mut source_buckets = FxHashMap::<BucketLoc, BucketData>::default();
            let mut target_buckets = FxHashMap::<BucketLoc, BucketData>::default();

            for observation in chunk {
                reos.insert(observation.reo_id);
                update_buckets(
                    observation,
                    &mut source_buckets,
                    &mut target_buckets,
                    &feature_buckets,
                );
            }

            (reos, source_buckets, target_buckets)
        })
        .collect();

    // Merge bucket collections together
    let mut reos = RoaringTreemap::new();
    let mut source_buckets = FxHashMap::<BucketLoc, BucketData>::default();
    let mut target_buckets = FxHashMap::<BucketLoc, BucketData>::default();

    for (rc, sb, tb) in filter_results {
        reos.extend(rc);
        update_bucket_map(&mut source_buckets, &sb);
        update_bucket_map(&mut target_buckets, &tb);
    }

    //
    // Build Final output data
    //

    // Turn bucket lists into chromosome
    let mut chromosomes: Vec<_> = if let Some(chromo_idx) = filters.chrom {
        data.chromosomes
            .iter()
            .filter(|c| c.index == chromo_idx)
            .map(|c| FilteredChromosome {
                chrom: c.chrom.clone(),
                index: c.index,
                bucket_size: data.bucket_size,
                target_intervals: Vec::new(),
                source_intervals: Vec::new(),
            })
            .collect()
    } else {
        data.chromosomes
            .iter()
            .map(|c| FilteredChromosome {
                chrom: c.chrom.clone(),
                index: c.index,
                bucket_size: data.bucket_size,
                target_intervals: Vec::new(),
                source_intervals: Vec::new(),
            })
            .collect()
    };

    let mut sources = RoaringTreemap::default();
    let mut targets = RoaringTreemap::default();
    let mut min_effect = f32::INFINITY;
    let mut max_effect = f32::NEG_INFINITY;

    let mut min_sig = f64::INFINITY;
    let mut max_sig = f64::NEG_INFINITY;

    gen_filtered_data(
        source_buckets,
        filters.chrom,
        &mut sources,
        &mut min_effect,
        &mut max_effect,
        &mut min_sig,
        &mut max_sig,
        &mut chromosomes
            .iter_mut()
            .map(|c| &mut c.source_intervals)
            .collect(),
        bucket_size,
        feature_buckets,
    );
    gen_filtered_data(
        target_buckets,
        filters.chrom,
        &mut targets,
        &mut min_effect,
        &mut max_effect,
        &mut min_sig,
        &mut max_sig,
        &mut chromosomes
            .iter_mut()
            .map(|c| &mut c.target_intervals)
            .collect(),
        bucket_size,
        feature_buckets,
    );

    // Make sure no numeric intervals include infinity
    min_effect = if min_effect == f32::INFINITY {
        effect_size_interval.0
    } else {
        min_effect
    };
    max_effect = if max_effect == f32::NEG_INFINITY {
        effect_size_interval.1
    } else {
        max_effect
    };

    min_sig = if min_sig == f64::INFINITY {
        sig_interval.0
    } else {
        min_sig
    };
    max_sig = if max_sig == f64::NEG_INFINITY {
        sig_interval.1
    } else {
        max_sig
    };

    FilteredData {
        chromosomes,
        bucket_size,
        numeric_intervals: FilterIntervals {
            effect: (min_effect, max_effect),
            sig: (min_sig, max_sig),
        },
        reo_count: reos.len(),
        sources,
        targets,
    }
}
