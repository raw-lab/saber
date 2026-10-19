//! Streaming FASTA reader: yields one record at a time with bounded memory.
//!
//! Unlike [`super::fasta::Fasta::from_file`], which materialises the whole
//! file, this reader keeps only the current record. It re-uses its buffers,
//! accepts multi-line sequences, CRLF line endings, blank lines and
//! lower-case residues, and validates what it can (sequence data before the
//! first header, empty identifiers).

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// One FASTA record borrowed from the reader's internal buffers.
pub struct FastaItem<'a> {
    /// First whitespace-delimited token of the header.
    pub id: &'a str,
    /// Full header line without the leading '>'.
    pub header: &'a str,
    /// Raw residues as found in the file (line breaks removed).
    pub seq: &'a [u8],
}

/// Streaming FASTA reader.
pub struct FastaStream<R: Read> {
    reader: BufReader<R>,
    line: Vec<u8>,
    header: String,
    next_header: Option<String>,
    seq: Vec<u8>,
    records: u64,
    done: bool,
}

impl FastaStream<File> {
    pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let f = File::open(path.as_ref()).map_err(|e| {
            std::io::Error::new(e.kind(), format!("{}: {}", path.as_ref().display(), e))
        })?;
        Ok(Self::new(f))
    }
}

impl<R: Read> FastaStream<R> {
    pub fn new(r: R) -> Self {
        Self {
            reader: BufReader::with_capacity(1 << 20, r),
            line: Vec::with_capacity(256),
            header: String::new(),
            next_header: None,
            seq: Vec::with_capacity(4096),
            records: 0,
            done: false,
        }
    }

    /// Number of records returned so far.
    pub fn records_read(&self) -> u64 {
        self.records
    }

    /// Advance to the next record. Returns `Ok(None)` at end of input.
    pub fn next_record(&mut self) -> std::io::Result<Option<FastaItem<'_>>> {
        if self.done {
            return Ok(None);
        }
        // Find the header for this record.
        let header = if let Some(h) = self.next_header.take() {
            h
        } else {
            loop {
                self.line.clear();
                let n = self.reader.read_until(b'\n', &mut self.line)?;
                if n == 0 {
                    self.done = true;
                    return Ok(None);
                }
                let l = trim_eol(&self.line);
                if l.is_empty() {
                    continue;
                }
                if l[0] == b'>' {
                    break String::from_utf8_lossy(&l[1..]).into_owned();
                }
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "FASTA: sequence data before the first '>' header",
                ));
            }
        };
        self.seq.clear();
        loop {
            self.line.clear();
            let n = self.reader.read_until(b'\n', &mut self.line)?;
            if n == 0 {
                self.done = true;
                break;
            }
            let l = trim_eol(&self.line);
            if l.is_empty() {
                continue;
            }
            if l[0] == b'>' {
                self.next_header = Some(String::from_utf8_lossy(&l[1..]).into_owned());
                break;
            }
            self.seq.extend(l.iter().filter(|b| !b.is_ascii_whitespace()));
        }
        self.header = header;
        self.records += 1;
        let id = self.header.split_whitespace().next().unwrap_or("");
        if id.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("FASTA: record {} has an empty identifier", self.records),
            ));
        }
        Ok(Some(FastaItem { id, header: &self.header, seq: &self.seq }))
    }
}

fn trim_eol(l: &[u8]) -> &[u8] {
    let mut end = l.len();
    while end > 0 && (l[end - 1] == b'\n' || l[end - 1] == b'\r') {
        end -= 1;
    }
    &l[..end]
}

/// Owned record, used when a batch of queries must outlive the reader.
#[derive(Debug, Clone)]
pub struct OwnedRecord {
    pub id: String,
    pub header: String,
    pub seq: Vec<u8>,
    /// 0-based ordinal in the input file (output ordering key).
    pub ordinal: u64,
}

impl<R: Read> FastaStream<R> {
    /// Read up to `max_records` records or `max_residues` residues (whichever
    /// comes first, but always at least one record if any remain).
    pub fn next_batch(&mut self, max_records: usize, max_residues: usize) -> std::io::Result<Vec<OwnedRecord>> {
        let mut out = Vec::new();
        let mut residues = 0usize;
        while out.len() < max_records && (out.is_empty() || residues < max_residues) {
            let ord = self.records;
            match self.next_record()? {
                Some(r) => {
                    residues += r.seq.len();
                    out.push(OwnedRecord {
                        id: r.id.to_string(),
                        header: r.header.to_string(),
                        seq: r.seq.to_vec(),
                        ordinal: ord,
                    });
                }
                None => break,
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiline_crlf_blank_lines() {
        let data = b"\n>a desc one\r\nMKV\r\nLAA\r\n\r\n>b\nGG\n>c  spaced id\n\nWW\nY";
        let mut s = FastaStream::new(&data[..]);
        let r = s.next_record().unwrap().unwrap();
        assert_eq!((r.id, r.header, r.seq), ("a", "a desc one", &b"MKVLAA"[..]));
        let r = s.next_record().unwrap().unwrap();
        assert_eq!((r.id, r.seq), ("b", &b"GG"[..]));
        let r = s.next_record().unwrap().unwrap();
        assert_eq!((r.id, r.seq), ("c", &b"WWY"[..]));
        assert!(s.next_record().unwrap().is_none());
        assert_eq!(s.records_read(), 3);
    }

    #[test]
    fn rejects_data_before_header_and_empty_ids() {
        let mut s = FastaStream::new(&b"MKV\n>a\nMM\n"[..]);
        assert!(s.next_record().is_err());
        let mut s = FastaStream::new(&b"> \nMM\n"[..]);
        assert!(s.next_record().is_err());
    }

    #[test]
    fn batches_respect_limits() {
        let data = b">a\nAAAA\n>b\nCCCC\n>c\nDDDD\n>d\nEE\n";
        let mut s = FastaStream::new(&data[..]);
        let b1 = s.next_batch(2, 1_000_000).unwrap();
        assert_eq!(b1.len(), 2);
        assert_eq!(b1[1].ordinal, 1);
        let b2 = s.next_batch(10, 5).unwrap();
        assert_eq!(b2.len(), 2); // stops once residues >= 5
        let b3 = s.next_batch(10, 5).unwrap();
        assert_eq!(b3.len(), 0);
    }
}
