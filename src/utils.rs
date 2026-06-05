/// Utility functions for progress tracking and statistics
use std::time::Instant;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// Progress bar for tracking work
pub struct ProgressBar {
    total: usize,
    current: Arc<AtomicUsize>,
    start: Instant,
    name: String,
}

impl ProgressBar {
    /// Create a new progress bar
    pub fn new(total: usize, name: &str) -> Self {
        Self {
            total,
            current: Arc::new(AtomicUsize::new(0)),
            start: Instant::now(),
            name: name.to_string(),
        }
    }

    /// Increment progress
    pub fn inc(&self) {
        self.current.fetch_add(1, Ordering::Relaxed);
    }

    /// Add to progress
    pub fn add(&self, n: usize) {
        self.current.fetch_add(n, Ordering::Relaxed);
    }

    /// Get current progress
    pub fn current(&self) -> usize {
        self.current.load(Ordering::Relaxed)
    }

    /// Get elapsed time
    pub fn elapsed(&self) -> std::time::Duration {
        self.start.elapsed()
    }

    /// Get rate (items per second)
    pub fn rate(&self) -> f64 {
        let elapsed_secs = self.start.elapsed().as_secs_f64();
        if elapsed_secs > 0.0 {
            self.current() as f64 / elapsed_secs
        } else {
            0.0
        }
    }

    /// Get eta in seconds
    pub fn eta(&self) -> f64 {
        let rate = self.rate();
        if rate > 0.0 {
            (self.total - self.current()) as f64 / rate
        } else {
            0.0
        }
    }

    /// Print progress
    pub fn print(&self) {
        let current = self.current();
        let percent = if self.total > 0 {
            (current as f64 / self.total as f64) * 100.0
        } else {
            0.0
        };
        let rate = self.rate();
        let eta = self.eta();

        eprintln!(
            "[{}] Progress: {}/{} ({:.1}%) Rate: {:.1} items/s ETA: {:.0}s",
            self.name, current, self.total, percent, rate, eta
        );
    }
}

/// Performance statistics
#[derive(Debug, Clone)]
pub struct Statistics {
    pub queries_processed: usize,
    pub total_alignments: usize,
    pub high_confidence: usize, // evalue < 0.001
    pub medium_confidence: usize, // 0.001 <= evalue < 0.1
    pub low_confidence: usize, // evalue >= 0.1
    pub mean_score: f64,
    pub median_score: f64,
    pub mean_identity: f64,
    pub elapsed_time: std::time::Duration,
}

impl Statistics {
    /// Create empty statistics
    pub fn new() -> Self {
        Self {
            queries_processed: 0,
            total_alignments: 0,
            high_confidence: 0,
            medium_confidence: 0,
            low_confidence: 0,
            mean_score: 0.0,
            median_score: 0.0,
            mean_identity: 0.0,
            elapsed_time: std::time::Duration::ZERO,
        }
    }

    /// Get throughput (alignments per second)
    pub fn throughput(&self) -> f64 {
        let secs = self.elapsed_time.as_secs_f64();
        if secs > 0.0 {
            self.total_alignments as f64 / secs
        } else {
            0.0
        }
    }

    /// Print statistics summary
    pub fn print_summary(&self) {
        println!("\n========== SABER Alignment Statistics ==========");
        println!("Queries processed:        {}", self.queries_processed);
        println!("Total alignments:         {}", self.total_alignments);
        println!("High confidence (evalue < 0.001): {}", self.high_confidence);
        println!("Medium confidence:        {}", self.medium_confidence);
        println!("Low confidence:           {}", self.low_confidence);
        println!("Mean alignment score:     {:.2}", self.mean_score);
        println!("Median alignment score:   {:.2}", self.median_score);
        println!("Mean identity:            {:.2}%", self.mean_identity);
        println!("Elapsed time:             {:.2}s", self.elapsed_time.as_secs_f64());
        println!("Throughput:               {:.0} alignments/s", self.throughput());
        println!("================================================\n");
    }
}

/// Timer utility for benchmarking
pub struct Timer {
    start: Option<Instant>,
    name: String,
}

impl Timer {
    /// Create a new timer
    pub fn new(name: &str) -> Self {
        Self {
            start: None,
            name: name.to_string(),
        }
    }

    /// Start timing
    pub fn start(&mut self) {
        self.start = Some(Instant::now());
    }

    /// Get elapsed time without stopping
    pub fn elapsed(&self) -> Option<std::time::Duration> {
        self.start.map(|s| s.elapsed())
    }

    /// Stop and print elapsed time
    pub fn stop_and_print(&self) {
        if let Some(elapsed) = self.elapsed() {
            eprintln!(
                "[{}] Completed in {:.3}s",
                self.name,
                elapsed.as_secs_f64()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_progress_bar() {
        let progress = ProgressBar::new(100, "test");
        for _ in 0..10 {
            progress.inc();
        }
        assert_eq!(progress.current(), 10);
    }

    #[test]
    fn test_statistics() {
        let stats = Statistics::new();
        assert_eq!(stats.queries_processed, 0);
        assert_eq!(stats.total_alignments, 0);
    }

    #[test]
    fn test_timer() {
        let mut timer = Timer::new("test");
        timer.start();
        thread::sleep(Duration::from_millis(10));
        let elapsed = timer.elapsed().unwrap();
        assert!(elapsed.as_millis() >= 10);
    }
}
