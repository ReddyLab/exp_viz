use roaring::RoaringTreemap;
use rustc_hash::FxHashSet;

use crate::filter_data_structures::*;

fn merge_chromosomes(
    result_data: &Vec<FilteredData>,
    chromosomes: Vec<String>,
) -> Vec<FilteredChromosome> {
    if result_data.len() == 0 {
        return Vec::new();
    } else if result_data.len() == 1 {
        return result_data[0].chromosomes.clone();
    }

    let mut new_coverage = Vec::new();
    let mut chroms_covered: FxHashSet<String> = FxHashSet::default();

    for chrom in chromosomes {
        // The default new_chromosome should never ever be used.
        let mut new_chromosome = FilteredChromosome {
            chrom: chrom.clone(),
            index: 0,
            bucket_size: 0,
            target_intervals: Vec::new(),
            source_intervals: Vec::new(),
        };
        for filtered_data in result_data {
            for filtered_chrom in &filtered_data.chromosomes {
                // Wrong chromosome
                if chrom != filtered_chrom.chrom {
                    continue;
                }

                // Right chromosome, but hasn't been added yet so no merging needed
                if !chroms_covered.contains(&filtered_chrom.chrom) {
                    chroms_covered.insert(filtered_chrom.chrom.clone());
                    new_chromosome = filtered_chrom.clone();
                    break;
                }

                // Right chromosome; has been added already; need to merge
                let mut source_intervals: Vec<FilteredBucket> = Vec::new();
                let mut i = 0; // new_chromosome source interval index
                let mut j = 0; // filtered_chrom source interval inded
                loop {
                    // We've run out of new_chromosome source intervals
                    // So we can just add the rest of the filtered_chrom source intervals
                    if i >= new_chromosome.source_intervals.len() {
                        while j < filtered_chrom.source_intervals.len() {
                            source_intervals.push(filtered_chrom.source_intervals[j].clone());
                            j += 1;
                        }
                        break;
                    }

                    // We've run out of filtered_chrom source intervals
                    // So we can just add the rest of the new_chromosome source intervals
                    if j >= filtered_chrom.source_intervals.len() {
                        while i < new_chromosome.source_intervals.len() {
                            source_intervals.push(new_chromosome.source_intervals[i].clone());
                            i += 1;
                        }
                        break;
                    }

                    if new_chromosome.source_intervals[i].start
                        < filtered_chrom.source_intervals[j].start
                    {
                        source_intervals.push(new_chromosome.source_intervals[i].clone());
                        i += 1;
                    } else if new_chromosome.source_intervals[i].start
                        > filtered_chrom.source_intervals[j].start
                    {
                        source_intervals.push(filtered_chrom.source_intervals[j].clone());
                        j += 1;
                    } else {
                        let mut assoc_buckets = filtered_chrom.source_intervals[j]
                            .associated_buckets
                            .clone();
                        assoc_buckets
                            .extend(new_chromosome.source_intervals[i].associated_buckets.iter());
                        source_intervals.push(FilteredBucket {
                            start: filtered_chrom.source_intervals[j].start,
                            count: filtered_chrom.source_intervals[j].count
                                + new_chromosome.source_intervals[i].count,
                            associated_buckets: assoc_buckets,
                            max_log10_sig: filtered_chrom.source_intervals[j]
                                .max_log10_sig
                                .max(new_chromosome.source_intervals[i].max_log10_sig),
                            max_abs_effect: if filtered_chrom.source_intervals[j]
                                .max_abs_effect
                                .abs()
                                > new_chromosome.source_intervals[i].max_abs_effect.abs()
                            {
                                filtered_chrom.source_intervals[j].max_abs_effect
                            } else {
                                new_chromosome.source_intervals[i].max_abs_effect
                            },
                        });
                        i += 1;
                        j += 1;
                    }
                }

                let mut target_intervals: Vec<FilteredBucket> = Vec::new();
                i = 0; // new_chromosome target interval index
                j = 0; // filtered_chrom target interval inded
                loop {
                    // We've run out of new_chromosome target intervals
                    // So we can just add the rest of the filtered_chrom target intervals
                    if i >= new_chromosome.target_intervals.len() {
                        while j < filtered_chrom.target_intervals.len() {
                            target_intervals.push(filtered_chrom.target_intervals[j].clone());
                            j += 1;
                        }
                        break;
                    }

                    // We've run out of filtered_chrom target intervals
                    // So we can just add the rest of the new_chromosome target intervals
                    if j >= filtered_chrom.target_intervals.len() {
                        while i < new_chromosome.target_intervals.len() {
                            target_intervals.push(new_chromosome.target_intervals[i].clone());
                            i += 1;
                        }
                        break;
                    }

                    if new_chromosome.target_intervals[i].start
                        < filtered_chrom.target_intervals[j].start
                    {
                        target_intervals.push(new_chromosome.target_intervals[i].clone());
                        i += 1;
                    } else if new_chromosome.target_intervals[i].start
                        > filtered_chrom.target_intervals[j].start
                    {
                        target_intervals.push(filtered_chrom.target_intervals[j].clone());
                        j += 1;
                    } else {
                        let mut assoc_buckets = filtered_chrom.target_intervals[j]
                            .associated_buckets
                            .clone();
                        assoc_buckets
                            .extend(new_chromosome.target_intervals[i].associated_buckets.iter());
                        target_intervals.push(FilteredBucket {
                            start: filtered_chrom.target_intervals[j].start,
                            count: filtered_chrom.target_intervals[j].count
                                + new_chromosome.target_intervals[i].count,
                            associated_buckets: assoc_buckets,
                            max_log10_sig: filtered_chrom.target_intervals[j]
                                .max_log10_sig
                                .max(new_chromosome.target_intervals[i].max_log10_sig),
                            max_abs_effect: if filtered_chrom.target_intervals[j]
                                .max_abs_effect
                                .abs()
                                > new_chromosome.target_intervals[i].max_abs_effect.abs()
                            {
                                filtered_chrom.target_intervals[j].max_abs_effect
                            } else {
                                new_chromosome.target_intervals[i].max_abs_effect
                            },
                        });
                        i += 1;
                        j += 1;
                    }
                }

                new_chromosome.source_intervals = source_intervals;
                new_chromosome.target_intervals = target_intervals;

                break;
            }
        }
        new_coverage.push(new_chromosome);
    }

    new_coverage
}

pub fn merge_filtered_data(
    result_data: Vec<FilteredData>,
    chromosome_list: Vec<String>,
) -> FilteredData {
    let chromosomes: Vec<FilteredChromosome> = merge_chromosomes(&result_data, chromosome_list);
    let numeric_intervals = result_data.iter().map(|d| d.numeric_intervals).fold(
        FilterIntervals {
            effect: (f32::MAX, f32::MIN),
            sig: (f64::MAX, f64::MIN),
        },
        |acc, d| FilterIntervals {
            effect: (acc.effect.0.min(d.effect.0), acc.effect.1.max(d.effect.1)),
            sig: (acc.sig.0.min(d.sig.0), acc.sig.1.max(d.sig.1)),
        },
    );

    FilteredData {
        chromosomes,
        numeric_intervals,
        bucket_size: result_data[0].bucket_size,
        reo_count: result_data.iter().map(|f| f.reo_count).sum(),
        sources: result_data
            .iter()
            .fold(RoaringTreemap::default(), |mut acc, f| {
                acc.extend(&f.sources);
                acc
            }),
        targets: result_data
            .iter()
            .fold(RoaringTreemap::default(), |mut acc, f| {
                acc.extend(&f.targets);
                acc
            }),
    }
}
