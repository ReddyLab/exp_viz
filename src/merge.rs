use pyo3::prelude::*;
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

                new_chromosome.index = filtered_chrom.index;

                // Right chromosome; has been added already; need to merge
                let mut source_intervals: Vec<FilteredBucket> = Vec::new();
                let mut i = 0;
                let mut j = 0;
                loop {
                    if i >= new_chromosome.source_intervals.len() {
                        while j < filtered_chrom.source_intervals.len() {
                            source_intervals.push(filtered_chrom.source_intervals[j].clone());
                            j += 1;
                        }
                        break;
                    }

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
                        });
                        i += 1;
                        j += 1;
                    }
                }
                new_chromosome.source_intervals = source_intervals;

                let mut target_intervals: Vec<FilteredBucket> = Vec::new();
                i = 0;
                j = 0;
                loop {
                    if i >= new_chromosome.target_intervals.len() {
                        while j < filtered_chrom.target_intervals.len() {
                            target_intervals.push(filtered_chrom.target_intervals[j].clone());
                            j += 1;
                        }
                        break;
                    }

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
                        });
                        i += 1;
                        j += 1;
                    }
                }
                new_chromosome.target_intervals = target_intervals;
                break;
            }
        }
        new_coverage.push(new_chromosome);
    }

    new_coverage
}

#[pyfunction]
pub fn merge_filtered_data(
    result_data: Vec<FilteredData>,
    chromosome_list: Vec<String>,
) -> PyResult<FilteredData> {
    let chromosomes: Vec<FilteredChromosome> = merge_chromosomes(&result_data, chromosome_list);
    let continuous_intervals = result_data.iter().map(|d| d.continuous_intervals).fold(
        FilterIntervals {
            effect: (f32::MAX, f32::MIN),
            sig: (f32::MAX, f32::MIN),
        },
        |acc, d| FilterIntervals {
            effect: (
                f32::min(acc.effect.0, d.effect.0),
                f32::max(acc.effect.1, d.effect.1),
            ),
            sig: (f32::min(acc.sig.0, d.sig.0), f32::max(acc.sig.1, d.sig.1)),
        },
    );

    Ok(FilteredData {
        chromosomes,
        continuous_intervals,
        item_counts: result_data
            .get(0)
            .and_then(|fd| Some(fd.item_counts))
            .unwrap_or([1, 2, 3]),
    })
}
