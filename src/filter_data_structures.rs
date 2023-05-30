use rustc_hash::{FxHashMap, FxHashSet};

use cov_viz_ds::{BucketLoc, ChromosomeData, CoverageData, DbID};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

/// Wraps the coverage data type so it can be passed to Python
#[pyclass(name = "CoverageData")]
pub struct PyCoverageData {
    pub wraps: CoverageData,
}

#[pyclass]
#[derive(Debug)]
pub struct Filter {
    #[pyo3(get, set)]
    pub discrete_facets: FxHashSet<DbID>,
    #[pyo3(get, set)]
    pub continuous_intervals: Option<FilterIntervals>,
}

#[pymethods]
impl Filter {
    #[new]
    pub fn new() -> Self {
        Filter {
            discrete_facets: FxHashSet::default(),
            continuous_intervals: None,
        }
    }

    pub fn __str__(&self) -> String {
        format!("Discrete Effects: {:?}", self.discrete_facets)
    }
}

#[pyclass]
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct FilterIntervals {
    #[pyo3(get, set)]
    pub effect: (f32, f32),
    #[pyo3(get, set)]
    pub sig: (f32, f32),
}

#[pymethods]
impl FilterIntervals {
    #[new]
    pub fn new() -> Self {
        FilterIntervals {
            effect: (f32::NEG_INFINITY, f32::INFINITY),
            sig: (f32::NEG_INFINITY, f32::INFINITY),
        }
    }

    pub fn __str__(&self) -> String {
        format!(
            "Effect Size: {:?}, Significance: {:?}",
            self.effect, self.sig
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, FromPyObject)]
pub struct FilteredBucket {
    pub start: u32,
    pub count: usize,
    pub associated_buckets: Vec<u32>,
}

#[pyclass]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FilteredChromosome {
    pub chrom: String,
    pub index: u8,
    pub bucket_size: u32,
    pub target_intervals: Vec<FilteredBucket>,
    pub source_intervals: Vec<FilteredBucket>,
}

#[pyclass]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FilteredData {
    #[pyo3(get, set)]
    pub chromosomes: Vec<FilteredChromosome>,
    #[pyo3(get, set)]
    pub continuous_intervals: FilterIntervals,
    #[pyo3(get, set)]
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
            continuous_intervals: FilterIntervals::new(),
            item_counts: [0, 0, 0],
        }
    }
}

#[pymethods]
impl FilteredData {
    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
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
