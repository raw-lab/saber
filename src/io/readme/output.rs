//! Output writers for the alignment formats SABER supports.
//!
//! `M8Writer` is the workhorse — BLAST tabular, one alignment per row, the
//! 12 standard columns. `M0Writer` produces the human-readable BLAST blocks.
//! SAM/BAM/JSON/CSV are also here.

use serde::Serialize;
use std::io::{BufWriter, Result as IoResult, Write};

use crate::algorithm::AlignmentResult;

pub trait Writer: Send {
    fn write_alignment(&mut self, alignment: &AlignmentResult) -> IoResult<()>;
    fn write_header(&mut self) -> IoResult<()>;
    fn write_footer(&mut self) -> IoResult<()>;
    fn flush(&mut self) -> IoResult<()>;
}

/// Write a single tab row in the BLAST `-outfmt 6` / `-m 8` layout.
fn write_m8_row<W: Write>(w: &mut W, a: &AlignmentResult) -> IoResult<()> {
    writeln!(
        w,
        "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
        a.query_id,
        a.subject_id,
        a.identity,
        a.align_length,
        a.mismatches,
        a.gaps,
        a.query_start,
        a.query_end,
        a.subject_start,
        a.subject_end,
        format_evalue(a.evalue),
        a.bitscore,
    )
}

/// Format an E-value the way BLAST does: `2e-50`, `0.0`, etc.
fn format_evalue(e: f64) -> String {
    if !e.is_finite() {
        return "inf".to_string();
    }
    if e == 0.0 {
        return "0.0".to_string();
    }
    if e >= 0.001 && e < 10.0 {
        format!("{:.3}", e)
    } else {
        // scientific notation, one digit before the dot
        format!("{:.2e}", e)
    }
}

pub struct M8Writer<W: Write> {
    writer: BufWriter<W>,
}
impl<W: Write> M8Writer<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer) } }
}
impl<W: Write + Send> Writer for M8Writer<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> { write_m8_row(&mut self.writer, a) }
    fn write_header(&mut self) -> IoResult<()> { Ok(()) /* BLAST m8 has no header */ }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

pub struct M9Writer<W: Write> {
    writer: BufWriter<W>,
    wrote_columns_line: bool,
}
impl<W: Write> M9Writer<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer), wrote_columns_line: false } }
}
impl<W: Write + Send> Writer for M9Writer<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        if !self.wrote_columns_line {
            writeln!(self.writer, "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")?;
            self.wrote_columns_line = true;
        }
        write_m8_row(&mut self.writer, a)
    }
    fn write_header(&mut self) -> IoResult<()> {
        writeln!(self.writer, "# SABER {}", crate::VERSION)
    }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

/// BLAST `-m 0` block-style output.
pub struct M0Writer<W: Write> {
    writer: BufWriter<W>,
}
impl<W: Write> M0Writer<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer) } }
}
impl<W: Write + Send> Writer for M0Writer<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        let pct_gaps = if a.align_length > 0 { 100.0 * a.gaps as f64 / a.align_length as f64 } else { 0.0 };
        writeln!(self.writer, ">{}", a.subject_id)?;
        writeln!(self.writer, " Score = {:.1} bits ({}), Expect = {}", a.bitscore, a.score, format_evalue(a.evalue))?;
        writeln!(self.writer, " Identities = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
            a.matches, a.align_length.max(1), a.identity, a.gaps, a.align_length.max(1), pct_gaps)?;
        writeln!(self.writer)?;
        // Print blocks of 60 like BLAST.
        let block = 60usize;
        let qb = a.query_seq.as_bytes();
        let sb = a.subject_seq.as_bytes();
        let mb = a.match_string.as_bytes();
        let mut q = a.query_start;
        let mut s = a.subject_start;
        let mut i = 0usize;
        while i < qb.len() {
            let j = (i + block).min(qb.len());
            let q_slice = &qb[i..j];
            let s_slice = &sb[i..j];
            let m_slice = &mb[i..j];
            let q_consumed = q_slice.iter().filter(|c| **c != b'-').count();
            let s_consumed = s_slice.iter().filter(|c| **c != b'-').count();
            writeln!(self.writer, "Query  {:>6}  {}  {}", q, String::from_utf8_lossy(q_slice), q + q_consumed.saturating_sub(1))?;
            writeln!(self.writer, "               {}", String::from_utf8_lossy(m_slice))?;
            writeln!(self.writer, "Sbjct  {:>6}  {}  {}", s, String::from_utf8_lossy(s_slice), s + s_consumed.saturating_sub(1))?;
            writeln!(self.writer)?;
            q += q_consumed;
            s += s_consumed;
            i = j;
        }
        Ok(())
    }
    fn write_header(&mut self) -> IoResult<()> {
        writeln!(self.writer, "SABER {}  https://github.com/raw937/saber", crate::VERSION)?;
        writeln!(self.writer)
    }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

pub struct SAMWriter<W: Write> { writer: BufWriter<W> }
impl<W: Write> SAMWriter<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer) } }
}
impl<W: Write + Send> Writer for SAMWriter<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        let mapq = if a.evalue.is_finite() && a.evalue < 0.001 { 60 } else { 30 };
        let cigar = if a.cigar.is_empty() { String::from("*") } else { a.cigar.clone() };
        let seq = if a.query_seq.is_empty() { String::from("*") } else { a.query_seq.replace('-', "") };
        writeln!(
            self.writer,
            "{}\t0\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tAS:i:{}\tNM:i:{}\tEV:f:{}",
            a.query_id, a.subject_id, a.subject_start, mapq, cigar, seq, a.score, a.mismatches, format_evalue(a.evalue),
        )
    }
    fn write_header(&mut self) -> IoResult<()> {
        writeln!(self.writer, "@HD\tVN:1.6\tSO:unsorted")?;
        writeln!(self.writer, "@PG\tID:saber\tPN:SABER\tVN:{}", crate::VERSION)
    }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

pub struct BAMWriter<W: Write> {
    writer: BufWriter<W>,
}
impl<W: Write> BAMWriter<W> {
    pub fn new(writer: W, _index: bool) -> Self { Self { writer: BufWriter::new(writer) } }
}
impl<W: Write + Send> Writer for BAMWriter<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        // BAM-as-text: SAM rows. Use real BAM (htslib) for production binary output.
        let mapq = if a.evalue.is_finite() && a.evalue < 0.001 { 60 } else { 30 };
        let cigar = if a.cigar.is_empty() { String::from("*") } else { a.cigar.clone() };
        let seq = if a.query_seq.is_empty() { String::from("*") } else { a.query_seq.replace('-', "") };
        writeln!(
            self.writer,
            "{}\t0\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tAS:i:{}\tNM:i:{}\tEV:f:{}",
            a.query_id, a.subject_id, a.subject_start, mapq, cigar, seq, a.score, a.mismatches, format_evalue(a.evalue),
        )
    }
    fn write_header(&mut self) -> IoResult<()> {
        writeln!(self.writer, "@HD\tVN:1.6\tSO:unsorted")?;
        writeln!(self.writer, "@PG\tID:saber\tPN:SABER\tVN:{}", crate::VERSION)
    }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

#[derive(Serialize)]
struct JsonAlignment<'a> {
    query_id: &'a str,
    subject_id: &'a str,
    score: i32,
    evalue: f64,
    bitscore: f64,
    identity: f64,
    align_length: usize,
    matches: usize,
    mismatches: usize,
    gaps: usize,
    query_range: (usize, usize),
    subject_range: (usize, usize),
    cigar: &'a str,
}

pub struct JSONWriter<W: Write> {
    writer: BufWriter<W>,
    first: bool,
}
impl<W: Write> JSONWriter<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer), first: true } }
}
impl<W: Write + Send> Writer for JSONWriter<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        // Streaming JSON array: emit comma separators between objects.
        if self.first {
            self.first = false;
        } else {
            writeln!(self.writer, ",")?;
        }
        let j = JsonAlignment {
            query_id: &a.query_id,
            subject_id: &a.subject_id,
            score: a.score,
            evalue: a.evalue,
            bitscore: a.bitscore,
            identity: a.identity,
            align_length: a.align_length,
            matches: a.matches,
            mismatches: a.mismatches,
            gaps: a.gaps,
            query_range: (a.query_start, a.query_end),
            subject_range: (a.subject_start, a.subject_end),
            cigar: &a.cigar,
        };
        let s = serde_json::to_string(&j).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        write!(self.writer, "  {}", s)
    }
    fn write_header(&mut self) -> IoResult<()> { writeln!(self.writer, "[") }
    fn write_footer(&mut self) -> IoResult<()> { writeln!(self.writer, "\n]") }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

pub struct CSVWriter<W: Write> { writer: BufWriter<W> }
impl<W: Write> CSVWriter<W> {
    pub fn new(writer: W) -> Self { Self { writer: BufWriter::new(writer) } }
}
impl<W: Write + Send> Writer for CSVWriter<W> {
    fn write_alignment(&mut self, a: &AlignmentResult) -> IoResult<()> {
        writeln!(
            self.writer,
            "{},{},{:.2},{},{},{},{},{},{},{},{},{:.1}",
            csv_escape(&a.query_id),
            csv_escape(&a.subject_id),
            a.identity, a.align_length, a.mismatches, a.gaps,
            a.query_start, a.query_end, a.subject_start, a.subject_end,
            format_evalue(a.evalue), a.bitscore,
        )
    }
    fn write_header(&mut self) -> IoResult<()> {
        writeln!(self.writer, "query_id,subject_id,identity,align_length,mismatches,gaps,q_start,q_end,s_start,s_end,evalue,bitscore")
    }
    fn write_footer(&mut self) -> IoResult<()> { Ok(()) }
    fn flush(&mut self) -> IoResult<()> { self.writer.flush() }
}

fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::AlignmentResult;

    fn sample() -> AlignmentResult {
        AlignmentResult {
            query_id: "q1".into(), subject_id: "s1".into(), score: 100,
            evalue: 1e-10, bitscore: 50.0, identity: 95.5,
            query_coverage: 100.0, subject_coverage: 100.0,
            query_start: 1, query_end: 100, subject_start: 1, subject_end: 100,
            align_length: 100, matches: 95, mismatches: 5, gaps: 0,
            cigar: "100M".into(),
            query_seq: "ACGT".into(), subject_seq: "ACGT".into(), match_string: "||||".into(),
        }
    }

    #[test]
    fn m8_has_12_tab_columns() {
        let mut out = Vec::new();
        {
            let mut w = M8Writer::new(&mut out);
            w.write_alignment(&sample()).unwrap();
            w.flush().unwrap();
        }
        let s = String::from_utf8(out).unwrap();
        let cols: Vec<&str> = s.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 12, "row: {:?}", s);
    }

    #[test]
    fn format_evalue_scientific_for_tiny() {
        assert!(format_evalue(1e-50).contains('e'));
    }

    #[test]
    fn format_evalue_decimal_for_midrange() {
        assert_eq!(format_evalue(0.5), "0.500");
    }

    #[test]
    fn json_writer_emits_array() {
        let mut out = Vec::new();
        {
            let mut w = JSONWriter::new(&mut out);
            w.write_header().unwrap();
            w.write_alignment(&sample()).unwrap();
            w.write_footer().unwrap();
            w.flush().unwrap();
        }
        let s = String::from_utf8(out).unwrap();
        assert!(s.starts_with('['));
        assert!(s.trim_end().ends_with(']'));
        assert!(s.contains("\"query_id\":\"q1\""));
    }
}
