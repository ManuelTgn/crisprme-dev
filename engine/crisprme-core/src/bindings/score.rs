use crate::common::iupac::Iupac;

#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("crisprme-core/include/scores.cuh");
        unsafe fn scores(
            query: *const u8,
            strings: *const u8,
            result: *mut u8,
            qlen: i32,
            slen: i32,
            n: i32,
        );
    }
}

pub fn scores_with_arena(
    query: &[Iupac],
    strings: &[Iupac],
    result: &mut [u8],
    slen: usize,
    n: usize,
) {
    assert!(query.len() <= slen);
    unsafe {
        ffi::scores(
            query.as_ptr() as *const u8,
            strings.as_ptr() as *const u8,
            result.as_mut_ptr(),
            query.len() as i32,
            slen as i32,
            n as i32,
        );
    }
}
