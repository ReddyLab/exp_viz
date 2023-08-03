use std::iter::zip;
use std::time::Instant;

use rustc_hash::FxHashSet;

use crate::filter_data_structures::*;
use cov_viz_ds::{BucketLoc, ChromosomeData, CoverageData, DbID, FacetCoverage, FacetRange};

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

pub fn filter_coverage_data(filters: &Filter, data: &CoverageData) -> FilteredData {
    let now = Instant::now();

    let mut coverage_data_disc_facets: FxHashSet<DbID> = FxHashSet::default();
    for facet in data.facets.iter() {
        match &facet.values {
            Some(facets) => facets.keys().for_each(|key| {
                coverage_data_disc_facets.insert(*key);
            }),
            None => (),
        };
    }

    let coverage_data_disc_facets: FxHashSet<DbID> = coverage_data_disc_facets
        .intersection(&filters.categorical_facets)
        .cloned()
        .collect();

    let skip_disc_facet_check = coverage_data_disc_facets.is_empty();
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
        Some(c) => FacetRange(c.sig.0, c.sig.1),
        None => data
            .facets
            .iter()
            .find(|f| f.name == "Significance")
            .unwrap()
            .range
            .unwrap(),
    };

    let mut source_facets: Vec<FxHashSet<DbID>> = Vec::new();
    for facet in data.facets.iter().filter(|f| {
        f.facet_type == "FacetType.CATEGORICAL"
            && f.coverage
                .as_ref()
                .unwrap()
                .contains(&FacetCoverage::Source)
    }) {
        source_facets.push(FxHashSet::from_iter(
            facet.values.as_ref().unwrap().keys().map(|k| *k),
        ));
    }

    let mut target_facets: Vec<FxHashSet<DbID>> = Vec::new();
    for facet in data.facets.iter().filter(|f| {
        f.facet_type == "FacetType.CATEGORICAL"
            && f.coverage
                .as_ref()
                .unwrap()
                .contains(&FacetCoverage::Target)
    }) {
        target_facets.push(FxHashSet::from_iter(
            facet.values.as_ref().unwrap().keys().map(|k| *k),
        ));
    }

    let sf_with_selections: Vec<&FxHashSet<DbID>> = source_facets
        .iter()
        .filter(|f| !f.is_disjoint(&coverage_data_disc_facets))
        .collect();
    let tf_with_selections: Vec<&FxHashSet<DbID>> = target_facets
        .iter()
        .filter(|f| !f.is_disjoint(&coverage_data_disc_facets))
        .collect();

    let selected_sf: Vec<Vec<DbID>> = sf_with_selections
        .iter()
        .map(|f| (*f & &coverage_data_disc_facets).into_iter().collect())
        .collect();
    let selected_tf: Vec<Vec<DbID>> = tf_with_selections
        .iter()
        .map(|f| (*f & &coverage_data_disc_facets).into_iter().collect())
        .collect();

    let filtered_data: Vec<(
        FilteredChromosome,
        f32,
        f32,
        f32,
        f32,
        FxHashSet<DbID>,
        FxHashSet<DbID>,
        FxHashSet<DbID>,
    )> = data
        .chromosomes
        .iter()
        .map(
            |chromosome| -> (
                FilteredChromosome,
                f32,
                f32,
                f32,
                f32,
                FxHashSet<DbID>,
                FxHashSet<DbID>,
                FxHashSet<DbID>,
            ) {
                let mut reo_ids = FxHashSet::<DbID>::default();
                let mut source_ids = FxHashSet::<DbID>::default();
                let mut target_ids = FxHashSet::<DbID>::default();
                let mut min_effect = f32::INFINITY;
                let mut max_effect = f32::NEG_INFINITY;

                let mut min_sig = f32::INFINITY;
                let mut max_sig = f32::NEG_INFINITY;

                let mut chrom_data = FilteredChromosome {
                    chrom: chromosome.chrom.clone(),
                    index: chromosome.index,
                    bucket_size: chromosome.bucket_size,
                    target_intervals: Vec::new(),
                    source_intervals: Vec::new(),
                };
                let chrom_info: Vec<(&ChromosomeData, &usize)> =
                    zip(&data.chromosomes, &data.chrom_lengths).collect();
                let bucket_list = BucketList::new(&chrom_info, chromosome.bucket_size as usize);
                if skip_cont_facet_check && sf_with_selections.is_empty() {
                    // do no filtering
                    for interval in &chromosome.source_intervals {
                        let mut associated_bucket = FxHashSet::<BucketLoc>::default();
                        let mut max_interval_sig = f32::MIN; // the max -log10 significance
                        let mut max_interval_effect: f32 = 0.0; // the largest _absolute_ effect size
                        for feature in &interval.values {
                            source_ids.insert(feature.id);
                            associated_bucket.extend(feature.associated_buckets.clone());
                            for observation in &feature.facets {
                                reo_ids.insert(observation.reo_id);
                                max_interval_sig =
                                    max_interval_sig.max(-observation.significance.log10());
                                max_interval_effect =
                                    if max_interval_effect.abs() > observation.effect_size.abs() {
                                        max_interval_effect
                                    } else {
                                        observation.effect_size
                                    }
                            }
                        }
                        chrom_data.source_intervals.push(FilteredBucket {
                            start: interval.start,
                            count: interval.values.len(),
                            associated_buckets: associated_bucket.iter().fold(
                                Vec::new(),
                                |mut acc, b| {
                                    acc.push(b.chrom as u32);
                                    acc.push(b.idx);

                                    acc
                                },
                            ),
                            max_log10_sig: max_interval_sig,
                            max_abs_effect: max_interval_effect,
                        });
                    }
                } else {
                    for interval in &chromosome.source_intervals {
                        let sources = &interval.values;
                        let mut new_source_count: usize = 0;
                        let mut new_target_buckets = bucket_list.clone();

                        let mut max_interval_sig = f32::MIN; // the max -log10 significance
                        let mut max_interval_effect: f32 = 0.0; // the largest _absolute_ effect size

                        for source in sources {
                            let mut new_regeffects = false;
                            for observation in &source.facets {
                                if skip_disc_facet_check
                                    || selected_sf
                                        .iter()
                                        .all(|sf| !is_disjoint(&observation.facet_ids, sf))
                                {
                                    let obs_sig = -observation.significance.log10();
                                    min_effect = observation.effect_size.min(min_effect);
                                    max_effect = observation.effect_size.max(max_effect);
                                    min_sig = obs_sig.min(min_sig);
                                    max_sig = obs_sig.max(max_sig);

                                    if skip_cont_facet_check
                                        || (observation.effect_size >= effect_size_interval.0
                                            && observation.effect_size <= effect_size_interval.1
                                            && obs_sig >= sig_interval.0
                                            && obs_sig <= sig_interval.1)
                                    {
                                        if !new_regeffects {
                                            new_regeffects = true;
                                            source_ids.insert(source.id);
                                        }
                                        reo_ids.insert(observation.reo_id);
                                        max_interval_sig = max_interval_sig.max(obs_sig);
                                        max_interval_effect = if max_interval_effect.abs()
                                            > observation.effect_size.abs()
                                        {
                                            max_interval_effect
                                        } else {
                                            observation.effect_size
                                        };
                                    }
                                }
                            }
                            if new_regeffects {
                                new_target_buckets.insert_from(&source.associated_buckets);
                                new_source_count += 1;
                            }
                        }

                        if new_source_count > 0 {
                            chrom_data.source_intervals.push(FilteredBucket {
                                start: interval.start,
                                count: new_source_count,
                                associated_buckets: new_target_buckets.flat_list(),
                                max_log10_sig: max_interval_sig,
                                max_abs_effect: max_interval_effect,
                            });
                        }
                    }
                }

                if skip_cont_facet_check && tf_with_selections.is_empty() {
                    // do no filtering
                    for interval in &chromosome.target_intervals {
                        let mut associated_bucket = FxHashSet::<BucketLoc>::default();
                        let mut max_interval_sig = f32::MIN; // the max -log10 significance
                        let mut max_interval_effect: f32 = 0.0; // the largest _absolute_ effect size
                        for feature in &interval.values {
                            target_ids.insert(feature.id);
                            associated_bucket.extend(feature.associated_buckets.clone());
                            for observation in &feature.facets {
                                reo_ids.insert(observation.reo_id);
                                max_interval_sig =
                                    max_interval_sig.max(-observation.significance.log10());
                                max_interval_effect =
                                    if max_interval_effect.abs() > observation.effect_size.abs() {
                                        max_interval_effect
                                    } else {
                                        observation.effect_size
                                    };
                            }
                        }
                        chrom_data.target_intervals.push(FilteredBucket {
                            start: interval.start,
                            count: interval.values.len(),
                            associated_buckets: associated_bucket.iter().fold(
                                Vec::new(),
                                |mut acc, b| {
                                    acc.push(b.chrom as u32);
                                    acc.push(b.idx);

                                    acc
                                },
                            ),
                            max_log10_sig: max_interval_sig,
                            max_abs_effect: max_interval_effect,
                        });
                    }
                } else {
                    for interval in &chromosome.target_intervals {
                        let targets = &interval.values;
                        let mut new_target_count: usize = 0;
                        let mut new_source_buckets = bucket_list.clone();

                        let mut max_interval_sig = f32::MIN; // the max -log10 significance
                        let mut max_interval_effect: f32 = 0.0; // the largest _absolute_ effect size

                        for target in targets {
                            let mut new_regeffects = false;
                            for observation in &target.facets {
                                if skip_disc_facet_check
                                    || selected_tf
                                        .iter()
                                        .all(|tf| !is_disjoint(&observation.facet_ids, tf))
                                {
                                    let obs_sig = -observation.significance.log10();
                                    min_effect = observation.effect_size.min(min_effect);
                                    max_effect = observation.effect_size.max(max_effect);
                                    min_sig = obs_sig.min(min_sig);
                                    max_sig = obs_sig.max(max_sig);

                                    if skip_cont_facet_check
                                        || (observation.effect_size >= effect_size_interval.0
                                            && observation.effect_size <= effect_size_interval.1
                                            && obs_sig >= sig_interval.0
                                            && obs_sig <= sig_interval.1)
                                    {
                                        if !new_regeffects {
                                            target_ids.insert(target.id);
                                            new_regeffects = true;
                                        }
                                        reo_ids.insert(observation.reo_id);
                                        max_interval_sig = max_interval_sig.max(obs_sig);
                                        max_interval_effect = if max_interval_effect.abs()
                                            > observation.effect_size.abs()
                                        {
                                            max_interval_effect
                                        } else {
                                            observation.effect_size
                                        };
                                    }
                                }
                            }
                            if new_regeffects {
                                new_source_buckets.insert_from(&target.associated_buckets);
                                new_target_count += 1;
                            }
                        }

                        if new_target_count > 0 {
                            chrom_data.target_intervals.push(FilteredBucket {
                                start: interval.start,
                                count: new_target_count,
                                associated_buckets: new_source_buckets.flat_list(),
                                max_log10_sig: max_interval_sig,
                                max_abs_effect: max_interval_effect,
                            })
                        }
                    }
                }

                return (
                    chrom_data, min_effect, max_effect, min_sig, max_sig, reo_ids, source_ids,
                    target_ids,
                );
            },
        )
        .collect();

    let mut min_effect = f32::INFINITY;
    let mut max_effect = f32::NEG_INFINITY;

    let mut min_sig = f32::INFINITY;
    let mut max_sig = f32::NEG_INFINITY;

    for x in filtered_data.iter() {
        min_effect = x.1.min(min_effect);
        max_effect = x.2.max(max_effect);
        min_sig = x.3.min(min_sig);
        max_sig = x.4.max(max_sig);
    }

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

    min_sig = if min_sig == f32::INFINITY {
        sig_interval.0
    } else {
        min_sig
    };
    max_sig = if max_sig == f32::NEG_INFINITY {
        sig_interval.1
    } else {
        max_sig
    };

    let item_sets = filtered_data.iter().fold(
        [
            FxHashSet::<DbID>::default(),
            FxHashSet::<DbID>::default(),
            FxHashSet::<DbID>::default(),
        ],
        |mut acc, data| {
            acc[0].extend(data.5.clone());
            acc[1].extend(data.6.clone());
            acc[2].extend(data.7.clone());
            acc
        },
    );

    let new_data = FilteredData {
        chromosomes: filtered_data.into_iter().map(|x| x.0).collect(),
        numeric_intervals: FilterIntervals {
            effect: (min_effect, max_effect),
            sig: (min_sig, max_sig),
        },
        item_counts: [
            item_sets[0].len() as u64,
            item_sets[1].len() as u64,
            item_sets[2].len() as u64,
        ],
    };

    println!("Time to filter data: {}ms", now.elapsed().as_millis());
    new_data
}
