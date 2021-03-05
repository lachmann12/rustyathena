pub fn fixedmi(var_a: &Vec<i32>, var_b: &Vec<i32>, bin_probs: &mut [[f32; 50]; 50], bin_number: usize) -> f32 {
    
    let bin_size = (var_a.len() / bin_number) as i32;
    
    for k in 0..bin_number {
        for v in 0..bin_number {
            bin_probs[k][v] = 0.0;
        }
    }
    
    for l in 0..var_a.len() {
        let i = var_a[l]/bin_size;
        let j = var_b[l]/bin_size;
        bin_probs[i as usize][j as usize] += 1.0;
    }

    let mut joint_entropy = 0.0;
    let vlen = var_a.len() as f32;
    for k in 0..bin_number {
        for v in 0..bin_number {
            let prob = bin_probs[k][v]/vlen;
            if prob != 0.0 {
                joint_entropy -= prob*prob.log2();
            }
        }
    }

    // H(X) + H(Y) - H(X,Y)
    // H(X) is log2(bin_number) for X uniformly distributed
    return 2.0*(bin_number as f32).log2() - joint_entropy;
}