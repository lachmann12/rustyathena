
pub fn calc_bspline(_rank: i32, _datalength: i32, _bin_number: usize) -> Vec<f32> {
    
    //println!("Rank: {} - Datalength: {} - Degree: {} - Bins: {}", _rank, _datalength, _degree, _binNumber);
    let _degree = 3;
    let bin_size = _datalength as f32/(_bin_number as f32);

    let nt = _bin_number + _degree + 1;
    
    let mut ti = vec![0f32; nt as usize];
    for i in 0..ti.len() {
        ti[(i as usize)] = (i as f32)*bin_size as f32;
    }
    //println!("{:?}", ti);

    let mut n = vec![0f32; (nt as usize)+1];
    let mut i1 = 0;
    
    // find bin of rank position
    while _rank as f32 >= ti[(i1 as usize)] {
        i1 = i1 + 1;
    }

    let i = i1 - 1;
    n[i] = 1.0;
    
    for m in 2.._degree+1 { // basis functions calculation
        let mut jb = i as i32 - m as i32 + 1;
        if jb < 0 {
            jb = 0;
        }
        
        for j in (jb as usize)..(i+1) {
            n[j] = n[j] * (_rank as f32 - ti[j]) / (ti[j + m - 1] - ti[j])
                    + n[j + 1] * (ti[j + m] - _rank as f32)
                    / (ti[j + m] - ti[j + 1]);
        }
    }

    let mut newn = vec![0f32; _bin_number as usize];

    for b in 0..newn.len()-2 {
        newn[b+1] = n[b];
    }
    
    if _rank as f32 <= (_degree as f32-1.0)*bin_size {
        newn[0] = 1.0 - newn[0] - n[0] - n[1];
    }
    
    if _rank as f32 >= ((_bin_number as f32) - _degree as f32 + 1.0)*bin_size {
        newn[_bin_number as usize -1] = 1.0 - newn[_bin_number as usize -2] - newn[_bin_number as usize -3];
    }

    return newn;
}