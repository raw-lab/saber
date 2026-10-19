//! SABER search engine: partitioned seed index, bounded seed processing,
//! the database-streaming scan schedule (v1.3), ungapped/gapped extension
//! and statistics shared by the CLI and the library API.

pub mod crc32;
pub(crate) mod hugebuf;
pub mod index;
pub mod engine;
pub mod driver;
pub mod presets;
pub mod api;
