use std::time::{Instant};


mod cliparameters;
mod readwrite;
mod normalize;
mod splinemi;
mod fixedbandwidthmi;

fn main() {

    // load user input parameters
    let matches = cliparameters::cli();

    // IO parameters
    let input_file = matches.value_of("input-file").unwrap();
    let output_file = matches.value_of("output-file").unwrap();
    let column_delimiter = matches.value_of("column-delimiter").unwrap_or("\t");
    let regulator_file = matches.value_of("regulators").unwrap_or("pairwise");

    // MI parameters
    let mi_method = matches.value_of("method").unwrap_or("spline");
    let number_bins = matches.value_of("bin-number").unwrap_or("12").parse::<usize>().unwrap();
    let thread_count = matches.value_of("thread-number").unwrap_or("4").parse::<usize>().unwrap();
    let sample_cutoff = matches.value_of("sample-cutoff").unwrap_or("1000000").parse::<usize>().unwrap();

    // DPI parameters 
    let dpi_mode = matches.is_present("verbose");
    
    // debugging
    let verbose = matches.is_present("verbose");
    let test_sample_count = matches.value_of("test-sample-count").unwrap_or("12").parse::<usize>().unwrap();
    // -------------------------------- end cli -------------------------------------------------


    // read data file
    let (gene_names, arr) = readwrite::read_from_file(&input_file, sample_cutoff, &column_delimiter);
    
    // read regulator / TF gene file
    let mut regulators = vec![];
    if regulator_file != "pairwise" {
        regulators = readwrite::read_regulators(&regulator_file, &column_delimiter, gene_names.clone());
        if verbose {
            println!("Regulators: {}", regulators.len());
        }
    }
    
    // initialize spline lookup table
    let mut spline_probabilities = vec![vec![0f32; number_bins]; arr.len()];
    for i in 0..arr.len() as i32 {
        spline_probabilities[i as usize] = splinemi::calc_bspline(i, arr.len() as i32, number_bins);
    }

    if verbose {
        println!("Samples: {}\nGenes: {}", arr.len(), arr[0].len());
        println!("Values: {} {} {}", arr[0][0], arr[0][1], arr[0][2]);
    }

    // normalize data by rank-rank transform. This is equivalent to quantile
    // transform of columns followed by rank transform of rows
    let now = Instant::now();
    let arr = normalize::rank_rank_transform(arr);
    if verbose {
        println!("Rank Time: {}ms", now.elapsed().as_millis()); 
    }

    // test fixed bandwidth

    if mi_method == "fixed" {
        let mut bin_probs_t = [[0f32; 50]; 50];

        let now = Instant::now();
        for h in 0..30000 {
            let fmi = fixedbandwidthmi::fixedmi(&arr[1], &arr[3], &mut bin_probs_t, number_bins);
        }
        println!("MI Time: {:.2}s", (now.elapsed().as_millis() as f32/1000.0));
    }

    let now = Instant::now();

    if verbose { 
        println!("Threads: {}",thread_count); 
    }

    if mi_method == "spline" {
        let (row_idx, col_names, mi_scores) = splinemi::calculate_spline_mi(arr, spline_probabilities, number_bins, thread_count, gene_names.clone(), regulators);
        if verbose {
            println!("MI Time: {:.2}s", (now.elapsed().as_millis() as f32/1000.0));
        }
    
        let now = Instant::now();
        if let Err(e) = readwrite::write_to_file(output_file, gene_names, row_idx, mi_scores) {
            eprintln!("{}", e)
        }
        if verbose {
            println!("Print rank time: {}ms", now.elapsed().as_millis());
        }
    } // end testing

}
