//! Large arrays backed by transparent huge pages where available.
//!
//! The scan schedule makes billions of random accesses into a few
//! multi-megabyte arrays (query windows, seed tables, database residues).
//! With 4 KiB pages most of those accesses also miss the TLB, and a page walk
//! costs as much as the access itself; 2 MiB pages cover the same arrays with
//! a handful of TLB entries. On Linux a buffer of at least [`MIN_BYTES`] is
//! an anonymous mapping aligned to 2 MiB and advised `MADV_HUGEPAGE` (used by
//! the kernel when transparent huge pages are set to `always` or `madvise`);
//! smaller buffers, other systems and failed mappings use an ordinary vector.
//! Contents and behaviour are the same either way.

#![allow(unsafe_code)]

/// Huge-page size and alignment.
const HUGE: usize = 2 << 20;
/// Smaller buffers use an ordinary vector (the second-level TLB covers a few
/// MB of 4 KiB pages; a huge page only pays off for larger arrays, and the
/// last one is partly unused).
pub const MIN_BYTES: usize = 4 << 20;

/// A fixed-length array of `Copy` values (see the module documentation).
pub struct HugeBuf<T: Copy> {
    ptr: *mut T,
    len: usize,
    /// (mapping start, mapping length) when backed by an anonymous mapping.
    map: Option<(usize, usize)>,
    vec: Vec<T>,
}

// SAFETY: HugeBuf owns its memory exclusively (like Vec<T>).
unsafe impl<T: Copy + Send> Send for HugeBuf<T> {}
// SAFETY: shared access is read-only through &[T] (like Vec<T>).
unsafe impl<T: Copy + Sync> Sync for HugeBuf<T> {}

impl<T: Copy> HugeBuf<T> {
    /// `len` copies of `v`.
    pub fn filled(len: usize, v: T) -> Self {
        let mut b = Self::alloc(len).unwrap_or_else(|| {
            let mut vec = vec![v; len];
            let ptr = vec.as_mut_ptr();
            HugeBuf { ptr, len, map: None, vec }
        });
        if b.map.is_some() {
            b.as_mut_slice().fill(v);
        }
        b
    }

    /// The contents of `v` (moved, or copied into huge pages when large).
    pub fn from_vec(mut v: Vec<T>) -> Self {
        let len = v.len();
        match Self::alloc(len) {
            Some(mut b) => {
                b.as_mut_slice().copy_from_slice(&v);
                b
            }
            None => {
                let ptr = v.as_mut_ptr();
                HugeBuf { ptr, len, map: None, vec: v }
            }
        }
    }

    /// True when backed by an advised huge-page mapping.
    #[allow(dead_code)]
    pub fn is_huge(&self) -> bool {
        self.map.is_some()
    }

    #[cfg(target_os = "linux")]
    fn alloc(len: usize) -> Option<Self> {
        let bytes = len.checked_mul(std::mem::size_of::<T>())?;
        if bytes < MIN_BYTES || std::mem::align_of::<T>() > HUGE || std::env::var_os("SABER_NO_HUGEPAGES").is_some() {
            return None;
        }
        let span = bytes.div_ceil(HUGE) * HUGE;
        let map_len = span + HUGE;
        // SAFETY: anonymous private mapping; checked for failure.
        let addr = unsafe {
            libc::mmap(std::ptr::null_mut(), map_len, libc::PROT_READ | libc::PROT_WRITE, libc::MAP_PRIVATE | libc::MAP_ANONYMOUS, -1, 0)
        };
        if addr == libc::MAP_FAILED {
            return None;
        }
        let start = addr as usize;
        let aligned = start.div_ceil(HUGE) * HUGE;
        // Huge pages for the whole 2 MiB blocks only: the partial last block
        // keeps 4 KiB pages, so no more memory is touched than the contents.
        // SAFETY: [aligned, aligned + whole) lies inside the mapping. The
        // advice is a hint; failure leaves ordinary pages.
        let whole = bytes / HUGE * HUGE;
        if whole > 0 {
            unsafe {
                libc::madvise(aligned as *mut libc::c_void, whole, libc::MADV_HUGEPAGE);
            }
        }
        Some(HugeBuf { ptr: aligned as *mut T, len, map: Some((start, map_len)), vec: Vec::new() })
    }

    #[cfg(not(target_os = "linux"))]
    fn alloc(_len: usize) -> Option<Self> {
        None
    }

    pub fn as_slice(&self) -> &[T] {
        // SAFETY: ptr points to len initialized values owned by self.
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }

    pub fn as_mut_slice(&mut self) -> &mut [T] {
        // SAFETY: as above, with exclusive access.
        unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
    }
}

impl<T: Copy> Drop for HugeBuf<T> {
    fn drop(&mut self) {
        #[cfg(target_os = "linux")]
        if let Some((start, map_len)) = self.map {
            // SAFETY: the mapping was created by `alloc` and is unmapped once.
            unsafe {
                libc::munmap(start as *mut libc::c_void, map_len);
            }
        }
        let _ = &self.vec;
    }
}

impl<T: Copy> std::ops::Deref for HugeBuf<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        self.as_slice()
    }
}

impl<T: Copy> std::ops::DerefMut for HugeBuf<T> {
    fn deref_mut(&mut self) -> &mut [T] {
        self.as_mut_slice()
    }
}

impl<T: Copy> Default for HugeBuf<T> {
    fn default() -> Self {
        Self::from_vec(Vec::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn large_and_small_buffers_hold_their_contents() {
        for len in [0usize, 10, 3 << 20, (5 << 20) + 7] {
            let mut b = HugeBuf::filled(len, 0xABu8);
            assert_eq!(b.len(), len);
            assert!(b.iter().all(|&x| x == 0xAB));
            if len > 0 {
                b[len - 1] = 1;
                assert_eq!(b[len - 1], 1);
            }
            let v: Vec<u32> = (0..len as u32 / 4).collect();
            let h = HugeBuf::from_vec(v.clone());
            assert_eq!(&h[..], &v[..]);
        }
    }
}
