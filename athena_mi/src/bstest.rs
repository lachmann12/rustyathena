
pub fn calcBSpline(_rank: i32, _datalength: i32, k: usize, binNumber: usize) -> Vec<f32> {
    
    println!("Rank: {} - Datalength: {} - Degree: {} - Bins: {}", _rank, _datalength, k, binNumber);

    let mut binsize = _datalength/(binNumber as i32);

    let mut ti = vec![0i32; (binNumber as usize)+4];
    for i in 1..binNumber+4 {
        ti[(i as usize)] = (i as i32)*binsize;
    }
    println!("{:?}", ti);

    let mut N = vec![0f32; (binNumber+4 as usize)+1];
    let mut i1 = 0;
    
    // find bin of rank position
    while _rank >= ti[(i1 as usize)] {
        i1 = i1 + 1;
    }
    
    //println!("{:?}", ti);

    let mut i = i1 - 1;
    N[i] = 1.0;
    
    for m in 2..k { // basis functions calculation
        let mut jb = i - m + 1;
        if jb < 0 {
            jb = 0;
        }
        
        for j in jb..(i+1) {
            N[j] = (N[j] * ((_rank - ti[j]) as f32) / ((ti[j + m - 1] - ti[j]) as f32)
                    + N[j + 1] * ((ti[j + m] - _rank) as f32)
                    / ((ti[j + m] - ti[j + 1]) as f32)).max(0.0);
        }
    }
    
    println!("{:?}", N);

    let mut newn = vec![0f32; (binNumber as usize)];
    
    for b in 0..(newn.len()-2) {
        newn[b+1] = N[b];
    }
    
    // edge case for summing to one on the left
    if _rank <= 2*binsize {
        newn[0] = 1.0-N[0]-N[1];
    }
    
    // edge case for summing to one on the right
    if _rank >= ((binNumber as i32)-2)*binsize{
        newn[(binNumber as usize)-1] = 1.0-N[(binNumber as usize)-3]-N[(binNumber as usize)-4];
    }
    
    return newn;
}