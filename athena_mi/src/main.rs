use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::{Instant};

use threadpool::ThreadPool;
use std::sync::mpsc::channel;

use std::error::Error;
use std::cmp;

use rand::Rng;
use frank::Ranking;

use std::sync::{Arc, RwLock};
use pbr::ProgressBar;

mod splinemi;

// CLI
pub fn cli<'a>() -> clap::ArgMatches<'a> {
    clap::App::new("Athena-MI")
      .version("0.1.0")
      .author("Alexander Lachmann <alexander.lachmann@gmail.com>")
      .about("Calculate probabilistic mutual information using B-splines")
      .arg(
        clap::Arg::with_name("input-file")
          .value_name("FILE")
          .short("i")
          .long("input-file")
          .takes_value(true)
          .required(true)
          .help("The gene count input file with genes as rows and samples as columns")
      )
      .arg(
        clap::Arg::with_name("output-file")
          .value_name("OUTFILE")
          .short("o")
          .long("output-file")
          .takes_value(true)
          .required(true)
          .help("The pairwise MI output")
      )
      .arg(
        clap::Arg::with_name("column-delimiter")
          .value_name("DELIMITER")
          .short("d")
          .long("column-delimiter")
          .takes_value(true)
          .required(false)
          .help("The delimiter separating each column, defaults to tab")
      )
      .arg(
        clap::Arg::with_name("thread-number")
          .value_name("THREADS")
          .short("t")
          .long("thread-number")
          .takes_value(true)
          .required(false)
          .help("Number of parallel threads, defaults to 4")
      )
      .arg(
        clap::Arg::with_name("test-sample-count")
          .value_name("scount")
          .short("s")
          .long("sampletest ")
          .takes_value(true)
          .required(false)
          .help("Number of samples for test, defaults to 4")
      )
      .arg(
        clap::Arg::with_name("bin-number")
          .value_name("BINS")
          .short("b")
          .long("bin-number")
          .takes_value(true)
          .required(false)
          .help("Number of bins.")
      )
      .arg(
        clap::Arg::with_name("sample-cutoff")
          .value_name("cutoff")
          .short("c")
          .long("sample-cutoff")
          .takes_value(true)
          .required(false)
          .help("Number of samples")
      )
      .get_matches()
}

fn main() {
    let matches = cli();
    let input_file = matches.value_of("input-file").unwrap();
    let output_file = matches.value_of("output-file").unwrap();
    let column_delimiter = matches.value_of("column-delimiter").unwrap_or("\t");
    let thread_count = matches.value_of("thread-number").unwrap_or("4").parse::<usize>().unwrap();
    let test_sample_count = matches.value_of("test-sample-count").unwrap_or("12").parse::<usize>().unwrap();
    let sample_cutoff = matches.value_of("sample-cutoff").unwrap_or("1000000").parse::<usize>().unwrap();

    let now = Instant::now();

    // open the file
    let f = BufReader::new(File::open(input_file).unwrap());
    let nrows = f.lines().count()-1;
    let mut f = BufReader::new(File::open(input_file).unwrap());
    // read the first line and extract the number from it
    let mut num_line = String::new();
    f.read_line(&mut num_line).unwrap();
    let sample_names: Vec<&str> = num_line.trim().split(column_delimiter).collect();
    let ncols = sample_names.len();
    
    // preallocate the array and read the data into it
    let mut arr = vec![vec![0f32; nrows]; cmp::min(sample_cutoff, ncols)];
    
    let mut gene_names: Vec<String> = vec!["".to_string(); nrows];
    let mut rng = rand::thread_rng();

    for (i, line) in f.lines().enumerate() {
        let lw = line.unwrap();
        for (j, number) in lw.split(column_delimiter).enumerate() {
            let entry = number.trim();
            if cmp::min(sample_cutoff, ncols) > j {
                if j == 0 {
                    gene_names[i] = String::from(entry);
                }
                else {
                    arr[j-1][i] = entry.parse().unwrap();
                    arr[j-1][i] = arr[j-1][i] + rng.gen::<f32>()/1000.0;
                }
            }
        }
    }
    
    println!("Samples: {}\nGenes: {}", arr.len(), arr[0].len());

    let number_bins = matches.value_of("bin-number").unwrap_or("12").parse::<usize>().unwrap();

    let mut spline_probabilities = vec![vec![0f32; number_bins]; arr.len()];
    for i in 0..arr.len() as i32 {
        spline_probabilities[i as usize] = splinemi::calc_bspline(i, arr.len() as i32, number_bins);
    }

    println!("Values: {} {} {}", arr[0][0], arr[0][1], arr[0][2]);
    println!("Time: {}ms", now.elapsed().as_millis());

    let now = Instant::now();

    let mut arr_rank = vec![vec![0f64; arr[0].len()]; arr.len()];
    for i in 0..arr.len() {
        arr_rank[i] = arr[i].rrank();
    }

    let mut arr = vec![vec![0f64; arr.len()]; arr[0].len()];
    for i in 0..arr_rank.len() {
        for j in 0..arr_rank[0].len() {
            arr[j][i] = arr_rank[i][j] + rng.gen::<f64>()/100000.0;
        }
    }

    let mut arr_rank = vec![vec![0i32; arr[0].len()]; arr.len()];
    for i in 0..arr.len() {
        let rr = arr[i].rrank();
        let mut rankround = vec![0i32; rr.len()];
        for j in 0..rr.len() {
            rankround[j] = rr[j] as i32;
        }
        arr_rank[i] = rankround;
    }
    
    println!("Rank Time: {}ms", now.elapsed().as_millis());

    let now = Instant::now();
    
    let mut tmi = 0.0;
    
    println!("Threads: {}",thread_count);
    let pool = ThreadPool::new(thread_count);
    let (tx, rx) = channel();

    let lock_arr = Arc::new(RwLock::new(arr_rank.clone()));
    let lock_spline = Arc::new(RwLock::new(spline_probabilities.clone()));
    let lock_gene_names = Arc::new(RwLock::new(gene_names.clone()));

    for tg in 0..arr_rank.len() as usize {
        let tx = tx.clone();
        let c_lock_arr = Arc::clone(&lock_arr);
        let c_lock_spline = Arc::clone(&lock_spline);
        let c_lock_gene_names = Arc::clone(&lock_gene_names);
        
        pool.execute(move || {
            let arr_rank = c_lock_arr.read().unwrap();
            let spline_prob = c_lock_spline.read().unwrap();
            let gene_names = c_lock_gene_names.read().unwrap();

            let mut bin_probs = [[0f32; 50]; 50];
            let mut mi = vec![0f32; gene_names.len()];

            let gene_a = &arr_rank[tg];
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
                tmi = mi_compute(&a_splines_all, &a_prob_all, &a_index_all, &arr_rank[g], &spline_prob, &mut bin_probs, number_bins as usize);
                mi[g] = tmi;
            }
            tx.send((tg, mi.clone())).expect("Could not send data!");
        });
    }

    let mut mi_results = vec![];
    let mut mi_genes = vec![];
    let mut pb = ProgressBar::new(arr_rank.len() as u64);

    drop(tx);
    for t in rx.iter() {
        pb.inc();
        let (tg, mi) = t;
        mi_genes.push(tg);
        mi_results.push(mi.clone());
    }
    pb.finish_print("Completed MI calculations");
    println!("MI Time: {:.2}s", (now.elapsed().as_millis() as f32/1000.0));

    let now = Instant::now();
    if let Err(e) = write_to_file(output_file, gene_names, mi_genes, mi_results) {
        eprintln!("{}", e)
    }
    println!("Print rank time: {}ms", now.elapsed().as_millis());
}

fn mi_compute(a_splines_all: &Vec<Vec<usize>>, a_prob_all: &Vec<Vec<f32>>, a_index_all: &Vec<usize>, gene_b: &Vec<i32>, spline_prob: &Vec<Vec<f32>>, bin_probs: &mut [[f32; 50]; 50], bin_number: usize) -> f32 {
    
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

fn retrieve_spline_probs(g1: i32, vv1: &mut [usize], vv2: &mut [f32], spline_prob: &Vec<Vec<f32>>) -> usize {
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

fn write_to_file(path: &str, genes: Vec<String>, mi_genes: Vec<usize>, data: Vec<Vec<f32>>) -> Result<(), Box<dyn Error>> {
    
    let f = File::create(path).expect("Unable to create file");
    let mut f = BufWriter::new(f);

    for i in 0..mi_genes.len() {
        if i > 0 {
            f.write_all("\t".as_bytes()).expect("Unable to write data");
        }
        f.write_all(genes[mi_genes[i]].as_bytes()).expect("Unable to write data");
    }
    f.write_all("\n".as_bytes()).expect("Unable to write data");

    for i in 0..genes.len() {
        f.write_all(genes[i].as_bytes()).expect("Unable to write data");
        for j in 0..mi_genes.len() {
            f.write_all("\t".as_bytes()).expect("Unable to write data");
            f.write_all((data[j][i].to_string()).as_bytes()).expect("Unable to write data");
        }
        f.write_all("\n".as_bytes()).expect("Unable to write data");
    }

    Ok(())
}