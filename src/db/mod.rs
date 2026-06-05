/// Database indexing and searching functionality
pub mod index;
pub mod search;

pub use index::{Database, DatabaseIndex};
pub use search::DatabaseSearch;

#[cfg(test)]
mod tests {
    #[test]
    fn test_module_structure() {
        // Ensure modules compile
    }
}
