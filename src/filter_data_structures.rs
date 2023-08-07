use rustc_hash::{FxHashMap, FxHashSet};

use cov_viz_ds::{BucketLoc, ChromosomeData, CoverageData, DbID};
use serde::{Deserialize, Serialize};

// When filtering this is the smallest we let a significance value be. Sometimes
// in they data the value is 0, which is infinity when we do a -log10 conversion,
// so we have to set an actual minimum. This number was selected as "resonable sounding".
// Don't be afraid to change it if another number becomes more "resonable sounding".
pub const MIN_SIG: f64 = 1e-100;

#[derive(Debug)]
pub struct Filter {
    pub categorical_facets: FxHashSet<DbID>,
    pub numeric_intervals: Option<FilterIntervals>,
}

impl Filter {
    pub fn new() -> Self {
        Filter {
            categorical_facets: FxHashSet::default(),
            numeric_intervals: None,
        }
    }

    pub fn __str__(&self) -> String {
        format!("Categorical Effects: {:?}", self.categorical_facets)
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct FilterIntervals {
    pub effect: (f32, f32),
    pub sig: (f64, f64),
}

impl FilterIntervals {
    pub fn new() -> Self {
        FilterIntervals {
            effect: (f32::NEG_INFINITY, f32::INFINITY),
            sig: (f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    pub fn __str__(&self) -> String {
        format!(
            "Effect Size: {:?}, Significance: {:?}",
            self.effect, self.sig
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FilteredBucket {
    pub start: u32,
    pub count: usize,
    pub associated_buckets: Vec<u32>,
    pub max_log10_sig: f64,  // Lower significance values are more significant
    pub max_abs_effect: f32, // largest absolute effect size
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FilteredChromosome {
    pub chrom: String,
    pub index: u8,
    pub bucket_size: u32,
    pub target_intervals: Vec<FilteredBucket>,
    pub source_intervals: Vec<FilteredBucket>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FilteredData {
    pub chromosomes: Vec<FilteredChromosome>,
    pub numeric_intervals: FilterIntervals,
    pub item_counts: [u64; 3],
}

impl FilteredData {
    pub fn from(data: &CoverageData) -> Self {
        FilteredData {
            chromosomes: data
                .chromosomes
                .iter()
                .map(|c| FilteredChromosome {
                    chrom: c.chrom.clone(),
                    index: c.index,
                    bucket_size: c.bucket_size,
                    target_intervals: Vec::new(),
                    source_intervals: Vec::new(),
                })
                .collect(),
            numeric_intervals: FilterIntervals::new(),
            item_counts: [0, 0, 0],
        }
    }
}

#[derive(Clone)]
pub struct BucketList {
    pub buckets: FxHashMap<u8, Vec<u32>>,
}

impl BucketList {
    pub fn new(chrom_info: &Vec<(&ChromosomeData, &usize)>, bucket_size: usize) -> Self {
        BucketList {
            buckets: chrom_info
                .iter()
                .map(|c| (c.0.index, vec![0; (c.1 / bucket_size) + 1]))
                .collect(),
        }
    }

    pub fn insert(&mut self, chrom: u8, bucket: usize) {
        let x = self.buckets.get_mut(&chrom);
        match x {
            Some(v) => v[bucket] = 1,
            None => (),
        };
    }

    pub fn insert_from<'a, I>(&mut self, from: I)
    where
        I: IntoIterator<Item = &'a BucketLoc>,
    {
        for bucket in from {
            self.insert(bucket.chrom, bucket.idx as usize);
        }
    }

    pub fn flat_list(&self) -> Vec<u32> {
        let mut new_list: Vec<u32> = Vec::new();
        for (i, chrom) in self.buckets.iter() {
            for (j, bucket) in chrom.iter().enumerate() {
                if *bucket == 1 {
                    new_list.push(*i as u32);
                    new_list.push(j as u32);
                }
            }
        }
        new_list
    }
}
