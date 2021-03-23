use std::cmp;
use std::fs::File;
use std::iter::FromIterator;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::error::Error;
use rand::Rng;

pub fn read_from_file(input_file: &str, sample_cutoff: usize, column_delimiter: &str) -> (Vec<String>, Vec<Vec<f32>>) {

    // open the file
    println!("File: {}", input_file);
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
                    arr[j-1][i] = entry.parse::<f32>().unwrap() + rng.gen::<f32>()/1000.0;
                }
            }
        }
    }

    return (gene_names, arr);
}


pub fn read_regulators(input_file: &str, column_delimiter: &str, gene_names: Vec<String>) -> Vec<String> {

    let f = BufReader::new(File::open(input_file).unwrap());
    
    let mut regulators: Vec<String> = vec![];
    
    for (i, line) in f.lines().enumerate() {
        let lw = line.unwrap();
        let v: Vec<&str> = lw.split(column_delimiter).collect();
        let reg = v[0].to_string();

        if gene_names.contains(&reg) {
            regulators.push(reg);
        }
    }
    
    return regulators;
}

pub fn write_to_file(path: &str, genes: Vec<String>, mi_genes: Vec<usize>, data: Vec<Vec<f32>>) -> Result<(), Box<dyn Error>> {
    
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