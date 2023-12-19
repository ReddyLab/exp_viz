use std::fmt;

use roaring::RoaringTreemap;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};

use cov_viz_ds::{BucketLoc, ChromosomeData, CoverageData, DbID};

// When filtering this is the smallest we let a significance value be. Sometimes
// in they data the value is 0, which is infinity when we do a -log10 conversion,
// so we have to set an actual minimum. This number was selected as "resonable sounding".
// Don't be afraid to change it if another number becomes more "resonable sounding".
pub const MIN_SIG: f64 = 1e-100;

#[derive(Debug)]
pub struct Filter {
    pub chrom: Option<u8>,
    pub categorical_facets: FxHashSet<DbID>,
    pub numeric_intervals: Option<FilterIntervals>,
}

impl Filter {
    pub fn new() -> Self {
        Filter {
            chrom: None,
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

#[derive(Clone, Debug)]
pub struct FilteredData {
    pub chromosomes: Vec<FilteredChromosome>,
    pub bucket_size: u32,
    pub numeric_intervals: FilterIntervals,
    pub reo_count: u64,
    pub sources: RoaringTreemap,
    pub targets: RoaringTreemap,
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
                    bucket_size: data.bucket_size,
                    target_intervals: Vec::new(),
                    source_intervals: Vec::new(),
                })
                .collect(),
            bucket_size: data.bucket_size,
            numeric_intervals: FilterIntervals::new(),
            reo_count: 0,
            sources: RoaringTreemap::default(),
            targets: RoaringTreemap::default(),
        }
    }
}

const FILTERED_DATA_CHROMOSOMES: &str = "chromosomes";
const FILTERED_DATA_BUCKET_SIZE: &str = "bucket_Size";
const FILTERED_DATA_NUMERIC_INTERVALS: &str = "numeric_intervals";
const FILTERED_DATA_REO_COUNT: &str = "reo_count";
const FILTERED_DATA_SOURCES: &str = "sources";
const FILTERED_DATA_TARGETS: &str = "targets";

impl Serialize for FilteredData {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("CoverageData", 2)?;
        state.serialize_field(FILTERED_DATA_CHROMOSOMES, &self.chromosomes)?;
        state.serialize_field(FILTERED_DATA_BUCKET_SIZE, &self.bucket_size)?;
        state.serialize_field(FILTERED_DATA_NUMERIC_INTERVALS, &self.numeric_intervals)?;
        state.serialize_field(FILTERED_DATA_REO_COUNT, &self.reo_count)?;
        let mut source_data = vec![];
        let _ = self.sources.serialize_into(&mut source_data);
        state.serialize_field(FILTERED_DATA_SOURCES, &source_data)?;
        let mut target_data = vec![];
        let _ = self.targets.serialize_into(&mut target_data);
        state.serialize_field(FILTERED_DATA_TARGETS, &target_data)?;

        state.end()
    }
}

impl<'de> Deserialize<'de> for FilteredData {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            Chromosomes,
            Bucket_Size,
            Numeric_Intervals,
            Reo_Count,
            Sources,
            Targets,
        }

        struct FilteredDataVisitor;

        impl<'de> Visitor<'de> for FilteredDataVisitor {
            type Value = FilteredData;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct FilteredData")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<FilteredData, V::Error>
            where
                V: SeqAccess<'de>,
            {
                let chromosomes = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let bucket_size = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let numeric_intervals = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                let reo_count = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let source_data: Vec<u8> = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                let target_data: Vec<u8> = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                let sources = RoaringTreemap::deserialize_from(&source_data[..]).unwrap();
                let targets = RoaringTreemap::deserialize_from(&target_data[..]).unwrap();

                Ok(FilteredData {
                    chromosomes,
                    bucket_size,
                    numeric_intervals,
                    reo_count,
                    sources,
                    targets,
                })
            }

            fn visit_map<V>(self, mut map: V) -> Result<FilteredData, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut chromosomes = None;
                let mut bucket_size = None;
                let mut numeric_intervals = None;
                let mut reo_count = None;
                let mut source_data: Option<Vec<u8>> = None;
                let mut target_data: Option<Vec<u8>> = None;
                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Chromosomes => {
                            if chromosomes.is_some() {
                                return Err(de::Error::duplicate_field(FILTERED_DATA_CHROMOSOMES));
                            }
                            chromosomes = Some(map.next_value()?);
                        }
                        Field::Bucket_Size => {
                            if bucket_size.is_some() {
                                return Err(de::Error::duplicate_field(FILTERED_DATA_BUCKET_SIZE));
                            }
                            bucket_size = Some(map.next_value()?);
                        }
                        Field::Numeric_Intervals => {
                            if numeric_intervals.is_some() {
                                return Err(de::Error::duplicate_field(
                                    FILTERED_DATA_NUMERIC_INTERVALS,
                                ));
                            }
                            numeric_intervals = Some(map.next_value()?);
                        }
                        Field::Reo_Count => {
                            if reo_count.is_some() {
                                return Err(de::Error::duplicate_field(FILTERED_DATA_REO_COUNT));
                            }
                            reo_count = Some(map.next_value()?);
                        }
                        Field::Sources => {
                            if source_data.is_some() {
                                return Err(de::Error::duplicate_field(FILTERED_DATA_SOURCES));
                            }
                            source_data = Some(map.next_value()?);
                        }
                        Field::Targets => {
                            if target_data.is_some() {
                                return Err(de::Error::duplicate_field(FILTERED_DATA_TARGETS));
                            }
                            target_data = Some(map.next_value()?);
                        }
                    }
                }
                let chromosomes = chromosomes
                    .ok_or_else(|| de::Error::missing_field(FILTERED_DATA_CHROMOSOMES))?;
                let bucket_size = bucket_size
                    .ok_or_else(|| de::Error::missing_field(FILTERED_DATA_BUCKET_SIZE))?;
                let numeric_intervals = numeric_intervals
                    .ok_or_else(|| de::Error::missing_field(FILTERED_DATA_NUMERIC_INTERVALS))?;
                let reo_count =
                    reo_count.ok_or_else(|| de::Error::missing_field(FILTERED_DATA_REO_COUNT))?;
                let source_data =
                    source_data.ok_or_else(|| de::Error::missing_field(FILTERED_DATA_SOURCES))?;
                let target_data =
                    target_data.ok_or_else(|| de::Error::missing_field(FILTERED_DATA_TARGETS))?;
                let sources = RoaringTreemap::deserialize_from(&source_data[..]).unwrap();
                let targets = RoaringTreemap::deserialize_from(&target_data[..]).unwrap();

                Ok(FilteredData {
                    chromosomes,
                    bucket_size,
                    numeric_intervals,
                    reo_count,
                    sources,
                    targets,
                })
            }
        }

        const FIELDS: &'static [&'static str] = &[
            FILTERED_DATA_CHROMOSOMES,
            FILTERED_DATA_BUCKET_SIZE,
            FILTERED_DATA_NUMERIC_INTERVALS,
            FILTERED_DATA_REO_COUNT,
            FILTERED_DATA_SOURCES,
            FILTERED_DATA_TARGETS,
        ];
        deserializer.deserialize_struct("FilteredData", FIELDS, FilteredDataVisitor)
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
