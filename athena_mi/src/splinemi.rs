use std::sync::{Arc, RwLock};
use threadpool::ThreadPool;
use std::sync::mpsc::channel;
use pbr::ProgressBar;

pub fn calc_bspline(_rank: i32, _datalength: i32, _bin_number: usize) -> Vec<f32> {
    
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


pub fn mi_compute(a_splines_all: &Vec<Vec<usize>>, a_prob_all: &Vec<Vec<f32>>, 
    a_index_all: &Vec<usize>, gene_b: &Vec<i32>, spline_prob: &Vec<Vec<f32>>, 
    bin_probs: &mut [[f32; 50]; 50], bin_number: usize) -> f32 {
    
    for k in 0..bin_number {
        for v in 0..bin_number {
            bin_probs[k][v] = 0.0;
        }
    }

    let mut b_splines = [0 as usize; 3];
    let mut b_prob = [0f32; 3];

    for l in 0..gene_b.len() {
        let b_index = retrieve_spline_probs(gene_b[l]-1, &mut b_splines, &mut b_prob, spline_prob);
        let a_index = &a_index_all[l];
        let a_splines = &a_splines_all[l];
        let a_prob = &a_prob_all[l];
        for i in 0..*a_index {
            for j in 0..b_index {
                bin_probs[a_splines[i]][b_splines[j]] += a_prob[i]*b_prob[j];
            }
        }
    }
    
    let mut mi = 0.0;
    let sq = (1.0/bin_number as f32)*(1.0/bin_number as f32);
    let norm = gene_b.len() as f32;
    for k in 0..bin_number {
        for v in 0..bin_number {
            if bin_probs[k][v] != 0.0 {
                let bp = bin_probs[k][v]/norm;
                mi += bp*((bp/sq).log2());
            }
        }
    }

    return mi;
}

pub fn retrieve_spline_probs(g1: i32, vv1: &mut [usize], vv2: &mut [f32], spline_prob: &Vec<Vec<f32>>) -> usize {
    let mut index = 0;
    let mut p = 0.0;
    let sp = &spline_prob[g1 as usize];
    for i in 0..sp.len() {
        p = sp[i];
        if sp[i] != 0.0 {
            vv1[index] = i;
            vv2[index] = p;
            index = index+1;
        }
        else if index > 0 {
            break;
        }
    }
    return index;
}

pub fn calculate_spline_mi(arr: Vec<Vec<i32>>, 
                            spline_probabilities: Vec<Vec<f32>>,
                            number_bins: usize,
                            thread_count: usize,
                            gene_names: Vec<String>, 
                            regulators: Vec<String>) -> (Vec<usize>, Vec<String>, Vec<Vec<f32>>) {

    let pool = ThreadPool::new(thread_count);
    let (tx, rx) = channel();

    let lock_arr = Arc::new(RwLock::new(arr.clone()));
    let lock_spline = Arc::new(RwLock::new(spline_probabilities.clone()));
    let lock_gene_names = Arc::new(RwLock::new(gene_names.clone()));
    
    for tg in 0..arr.len() as usize {
        let tx = tx.clone();
        let c_lock_arr = Arc::clone(&lock_arr);
        let c_lock_spline = Arc::clone(&lock_spline);
        let c_lock_gene_names = Arc::clone(&lock_gene_names);
        
        pool.execute(move || {
            let arr = c_lock_arr.read().unwrap();
            let spline_prob = c_lock_spline.read().unwrap();
            let gene_names = c_lock_gene_names.read().unwrap();

            let mut bin_probs = [[0f32; 50]; 50];
            let mut mi = vec![0f32; gene_names.len()];

            let gene_a = &arr[tg];
            let mut a_splines = vec![0 as usize; 3];
            let mut a_prob = vec![0f32; 3];
            let mut a_splines_all = vec![vec![0 as usize; 3]; gene_a.len()];
            let mut a_prob_all = vec![vec![0f32; 3]; gene_a.len()];
            let mut a_index_all = vec![0 as usize; gene_a.len()];
            
            for l in 0..gene_a.len() {
                let a_index = retrieve_spline_probs(gene_a[l]-1, &mut a_splines, &mut a_prob, &spline_prob);
                a_splines_all[l] = a_splines.to_vec();
                a_prob_all[l] = a_prob.to_vec();
                a_index_all[l] = a_index;
            }
            for g in 0..gene_names.len() {
                mi[g] = mi_compute(&a_splines_all, &a_prob_all, &a_index_all, 
                    &arr[g], &spline_prob, &mut bin_probs, number_bins as usize);
            }
            tx.send((tg, mi.clone())).expect("Could not send data!");
        });
    }

    let mut mi_results = vec![];
    let mut mi_genes = vec![];
    let mut pb = ProgressBar::new(arr.len() as u64);

    drop(tx);
    for t in rx.iter() {
        pb.inc();
        let (tg, mi) = t;
        mi_genes.push(tg);
        mi_results.push(mi.clone());
    }
    pb.finish_print("Completed MI calculations");
    
    return (mi_genes, gene_names, mi_results);
}