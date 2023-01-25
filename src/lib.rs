use pyo3::prelude::*;
use probability::distribution::{Discrete, Binomial};
use rand::distributions::{Weighted, WeightedChoice, Distribution, Poisson};
use rand::prelude::*;
use csv;
use std::io::{BufReader, BufRead, BufWriter, Write};
use std::fs;
use intervaltree::*;
use std::collections::HashMap;
use core::ops::Range;
use rulinalg::vector::Vector;

pub struct Read {
    chr: String,
    start: usize,
    end: usize,
    m: usize,
    t: usize,
}
pub struct Atlas {
    region_map: HashMap<String, IntervalTree<usize, usize>>,
    values: Vec<Vec<f64>>,
    pub cell_types: Vec<String>
}
impl Atlas {
    pub fn new(atlas: &str) -> Atlas {
        // Read atlas into memory
        let mut atlas_reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(atlas).expect("Could not open bed file");
        let mut region_desc_by_chr = HashMap::<String, Vec<(Range<usize>, usize)>>::new();

        //Store headers in cell types
        let header = atlas_reader.headers().unwrap();
        let cell_types: Vec<String> = header.iter().skip(3).map(|x| x.to_string()).collect();

        // this stores the atlas data for each interval
        let mut atlas_data = Vec::new();

        for r in atlas_reader.records() {
            let record = r.expect("Could not parse bed record");
            let chr: String  = record[0].to_string();
            let start: usize = record[1].parse().unwrap();
            let end: usize = record[2].parse().unwrap();
            let values: Vec<f64> = record.iter().skip(3).take(39).map(|x| x.parse().unwrap_or(0.0)).collect();

            let region_desc = region_desc_by_chr.entry(chr).or_insert( Vec::new() );
            region_desc.push( (start..end, atlas_data.len()) );
            atlas_data.push(values);
        }
        // build chr -> intervaltree map
        // the intervaltree allows us to look up an interval_idx for a given chromosome and position
        let mut interval_trees = HashMap::<String, IntervalTree<usize, usize>>::new();
        for (chr, region_desc) in region_desc_by_chr {
            interval_trees.insert(chr, region_desc.iter().cloned().collect());
        }
        Atlas {
            region_map: interval_trees,
            values: atlas_data,
            cell_types: cell_types
        }
    }
    pub fn query(&self, read: &Read) -> &Vec<f64> {
        //Query the region_map to get the methylation propensity
        let tree = self.region_map.get(&read.chr).unwrap();
        let mut query = tree.query(read.start..read.end);
        //println!("{}\t{}\t{}\n{:?}", read.chr,read.start,read.end,query);
        //tree.iter().take(20).for_each(|x| println!("{:?}",x));
        let atlas_idx = query.next().unwrap().value;
        &self.values[atlas_idx]
    }
}

#[pyclass]
pub struct MMSE {
    #[pyo3(get)]
    pub sigma: Vec<f64>,
    K: usize,
    N: usize,
    pub atlas: Atlas,
    methylome: Vec<Read>,
    p01: f64,
    p11: f64,
    #[pyo3(get)]
    pub assignments: Vec<usize>
}
impl MMSE {
    fn get_p(&self, xi: f64) -> f64 {
        self.p11*xi + self.p01*(1.0 - xi)
    }
    fn gamma_tilde(&self, k: usize, p: f64, read: &Read, sigma: &Vec<f64>) -> f64 {
        let bin = Binomial::new(read.t, p);
        sigma[k] * bin.mass(read.m)
    }
    pub fn gamma(&self, k: usize, atlas_row: &Vec<f64>, read: &Read) -> f64 {
        let sum : f64 = (0..self.K)
            .map(|k0| self.gamma_tilde(k0, self.get_p(atlas_row[k0]), read, &self.sigma)).sum();
        self.gamma_tilde(k, self.get_p(atlas_row[k]), read, &self.sigma) / sum
    }
    fn e_step(&self, sigma: &Vec<f64>) -> f64{
        // compute the log likelihood of the data
        (0..self.N)
            .map(|i| (self.atlas.query(&self.methylome[i]), &self.methylome[i]) )
            .map(|(atlas_row, read)| (0..self.K)
                 .map(|k| self.gamma_tilde(k, self.get_p(atlas_row[k]), read, sigma) )
                 .sum::<f64>().ln() )
            .sum()
    }
    fn print_gamma(&self, i: usize) -> (){
        // compute the log likelihood of the data
        let (atlas_row, read) = (self.atlas.query(&self.methylome[i]), &self.methylome[i]);
        (0..self.K).for_each(|k| print!("{}:{:.1}\t",k,self.gamma(k, atlas_row, read).ln() ));
        println!("gamma_sum: {}",
                 (0..self.K).map(|k| self.gamma(k, atlas_row, read).ln() ).sum::<f64>());
    }
    fn m_step(&mut self){
        // Compute MLE and update sigma
        self.sigma = (0..self.K)
            .map(|k| (0..self.N)
                 .map(|i| (self.atlas.query(&self.methylome[i]), &self.methylome[i]) )
                 .map(|(atlas_row, read)| self.gamma(k, atlas_row, read))
                 .sum::<f64>() / self.N as f64 )
            .collect();
    }
    fn noise_trim(&mut self, t: f64) {
        for k in 0..self.K {
            if self.sigma[k] < t { self.sigma[k] = 0.0 }
        }
        let sum: f64 = self.sigma.iter().sum();
        for k in 0..self.K {
            self.sigma[k] = self.sigma[k]/sum;
        }
    }
}


#[pymethods]
impl MMSE {
    #[new]
    pub fn new(methylome: &str, atlas: &str, init_sigma: Vec<f64>, p01: f64, p11: f64) -> PyResult<Self> {
        let atlas = Atlas::new(atlas);
        // Read methylome into memory
        // Store the list of intervals for the methylome instead of interval tree
        // Make a struct for reads
        let mut methylome_reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(methylome).expect("Could not open bed file");

        // get header from methylome reader
        let header = methylome_reader.headers().unwrap();
        // Map header to index
        let header_map = header.iter().enumerate().map(|(i, x)| (x.to_string(), i)).collect::<HashMap<String, usize>>();
        // get the index of the columns we need
        let chr_idx = header_map.get("chr").unwrap();
        let start_idx = header_map.get("start").unwrap();
        let end_idx = header_map.get("end").unwrap();
        let m_idx = header_map.get("modified_calls").unwrap();
        let t_idx = header_map.get("total_calls").unwrap();
        // this stores the methylome data for each interval
        let mut methylome_data = Vec::new();
        for r in methylome_reader.records() {
            let record = r.expect("Could not parse bed record");
            let read = Read {
                chr: record[*chr_idx].to_string(),
                start: record[*start_idx].parse().unwrap(),
                end: record[*end_idx].parse().unwrap(),
                t: record[*t_idx].parse().unwrap(),
                m: record[*m_idx].parse().unwrap(),
            };
            methylome_data.push(read);
        }
        let K = atlas.values[0].len();
        let N = methylome_data.len();

        Ok(MMSE{
            sigma: init_sigma,
            N: N,
            K: K,
            atlas: atlas,
            methylome: methylome_data,
            p01: p01,
            p11: p11,
            assignments: vec![0 as usize; N]})
    }
    pub fn optimize(&mut self, stop_thresh: f64, max_iter: u32, min_proportion: f64)-> PyResult<()>{
        let mut i = 0;
        let mut ll_prev = self.e_step(&self.sigma);
        let mut ll = self.e_step(&self.sigma);
        eprintln!("Running EM Inference");
        loop {
            i += 1;
            self.m_step();
            ll = self.e_step(&self.sigma);
            let condition = 100.0*((ll - ll_prev)/ll_prev).abs();
            eprintln!("Iteration: {}\tLog-likelihood: {:.5}\tPercent_change: {:.8}",
                      i, ll, condition);
            if ( condition < stop_thresh && i > 10 ) 
                || i == max_iter-1 {
                break;
            } else {
                ll_prev = ll;
            }
        }
        self.noise_trim(min_proportion);
        self.m_step();
        ll = self.e_step(&self.sigma);
        i += 1;
        eprintln!("Model optimized in {} steps. Log-likelihood = {}", i, ll);
        eprintln!("log-likelihood:\t{:.2}", ll);
        Ok(())
    }
    pub fn evaluate(&mut self, stop_thresh: f64, max_iter: u32, min_proportion: f64,
                    sigma: Vec<f64>, target_assignments: Vec<usize>) -> PyResult<()>{
        let mut i = 0;
        let mut ll_prev = self.e_step(&self.sigma);
        let mut ll = self.e_step(&self.sigma);
        loop {
            let loss = deconvolution_loss(&sigma, &self.sigma);
            let acc = self.accuracy(target_assignments.clone(), 0.5);
            eprintln!("Iter: {}\tLoss: {:.3}\tAccuracy: {:.3}\tLog-Likelihood: {:.3}\t True Log-likelihood: {:.3}",
                      i, loss, acc, ll, self.e_step(&sigma));
            i += 1;
            self.m_step();
            ll = self.e_step(&self.sigma);
            if (ll_prev - ll).abs()/ll < stop_thresh || i == max_iter {
                break;
            } else {
                ll_prev = ll;
            }
        }
        self.noise_trim(min_proportion);
        self.m_step();
        ll = self.e_step(&self.sigma);
        let loss = deconvolution_loss(&sigma, &self.sigma);
        let acc = self.accuracy(target_assignments.clone(), 0.5);
        eprintln!("Loss: {:.3}\t Accuracy: {:.3}\tLog-Likelihood: {:.3}\t True Log-likelihood: {:.3}", loss, acc, ll, self.e_step(&sigma));
        Ok(())
    }
    pub fn reset(&mut self){
        self.sigma = (0..self.K).map(|_k| 1.0/self.K as f64).collect()
    }
    fn assign_fragments(&mut self) -> () {
        for i in 0..self.N{ 
            let read = &self.methylome[i];
            let atlas_row = self.atlas.query(read);
            let gamma_vec: Vec<f64> = (0..self.K).map(|k| self.gamma(k, atlas_row, read)).collect();
            self.assignments[i] = Vector::new(gamma_vec).argmax().0;
        }
        ()
    }
    pub fn assign_fragments_t(&mut self, threshold: f64) -> Vec<usize> {
        // assign reads to the most likely cell type with confidence above threshold
        // If the read is unnassigned it gets value K
        let mut assignments = vec![0 as usize; self.N];
        for i in 0..self.N{ 
            let read = &self.methylome[i];
            let atlas_row = self.atlas.query(read);
            let gamma_vec: Vec<f64> = (0..self.K).map(|k| self.gamma(k, atlas_row, read)).collect();
            if gamma_vec.iter().cloned().fold(0./0., f64::max) > threshold {
                assignments[i] = Vector::new(gamma_vec).argmax().0;
            } else {
                // give the assignment a NaN value
                assignments[i] = self.K;
            }
        }
        assignments
    }
    pub fn accuracy(&mut self, target_assignments: Vec<usize>, threshold: f64) -> f64 {
        self.assign_fragments();
        self.assignments.iter().zip(target_assignments.iter())
            .map(|(value, target)| (*value == *target) as u32 as f64)
            .sum::<f64>() / self.assignments.len() as f64
    }
    pub fn cell_type_proportions(&mut self) -> HashMap<String, f64> {
        // Return a hashmap of keys = cell type names and values = respective sigma proportions
        let mut cell_type_proportions : HashMap<String, f64> = HashMap::new();
        for (k, cell_type) in self.atlas.cell_types.iter().enumerate() {
            cell_type_proportions.insert(cell_type.clone(), self.sigma[k]);
        }
        cell_type_proportions
    }
    pub fn class_probabilities(&mut self, target_assignments: Vec<usize>) -> Vec<Vec<f64>> {
        self.assign_fragments();
        let mut probs: Vec<Vec<f64>> = self.assignments.iter().zip(target_assignments.iter())
            .map(|(value, target)| (*value == *target) as u32 as f64 )
            .map(|correct| vec![correct] )
            .collect();
        (0..self.N).map(|i| (self.atlas.query(&self.methylome[i]), &self.methylome[i]) )
            .map(|(atlas_row, read)| (0..self.K).map(|k| self.gamma(k, atlas_row, read)).collect())
            .enumerate()
            .for_each(|(i, mut row)| probs[i].append(&mut row));
        (0..self.N)
            .map(|i| self.atlas.query(&self.methylome[i]).clone() ).enumerate()
            .for_each(|(i, mut row)| probs[i].append(&mut row) );
        probs
    }
}
pub fn deconvolution_loss(sigma: &Vec<f64>, sigma_hat: &Vec<f64>) -> f64 {
    sigma.iter().zip(sigma_hat.iter()).map(|(y, yhat)| (y-yhat).powf(2.0)).sum::<f64>().powf(0.5)
}
#[pyfunction]
pub fn generate_methylome(atlas: &str, sigma_path: &str, coverage: f64, region_size: u64, p01: f64, p11: f64) -> PyResult<()> {
    let atlas_file = fs::File::open(atlas)?;
    let mut atlas_reader = BufReader::new(atlas_file);
    let mut header = String::new();
    let _res = atlas_reader.read_line(&mut header);
    let cell_types: Vec<String> = header.split('\t')
        .filter(|x| !"chromosomestartendlabel".contains(x))
        .map(|x| x.trim())
        .map(|x| String::from(x))
        .collect();

    // open sigma tsv file and read it into a vector
    // warning: assumes that sigma.tsv lists the cell types in the same order as the atlas
    let contents = fs::read_to_string(sigma_path)?;
    let lines = contents.lines();
    let mut sigma: Vec<f64> = Vec::new();
    lines.skip(1).for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        sigma.push(fields[1].parse::<f64>().unwrap());
    });
    if cell_types.len() != sigma.len() { panic!("Length of Sigma ({}) != number of cell types ({})", sigma.len(), cell_types.len()) }

    let poi = Poisson::new(coverage);
    let mut rng = StdRng::from_entropy();
    let mut choices: Vec<Weighted<usize>> = sigma.iter().enumerate()
        .map(|(i, x)| Weighted { weight: (x*1000.0) as u32, item: i })
        .collect();
    let wc = WeightedChoice::new(&mut choices);

    println!("chromosome\tstart\tend\ttotal_calls\tmodified_calls\tcell_type\n");

    let lines: Vec<Vec<String>> = atlas_reader.lines()
            .filter_map(|line| line.ok())
            .map(|line| line.split('\t').map(|x| String::from(x)).collect())
            .collect();
    for line in lines {
        // Generate coverage from a poisson distribution
        let t = region_size; // rng.gen_range(5, 7);
        let c = poi.sample(&mut rng);
        for coverage in 0..c as usize {
            // Generate a random number of total calls and modified calls
            let cell_type = wc.sample(&mut rng);
            let z = &line[cell_type + 3].parse::<f64>().unwrap();
            //let p = z*p11 + (1.0-z)*p01;
            let bin = rand::distributions::Binomial::new(t, *z);
            let mut m = bin.sample(&mut rng);

            // Introduce errors at a rate of p01 bit switch on modified calls
            // and (1-p11) bit `switch on unmodified calls
            let errors_01 = rand::distributions::Binomial::new(t-m, p01).sample(&mut rng);
            let errors_10 = rand::distributions::Binomial::new(m, 1f64-p11).sample(&mut rng);
            m = m + errors_01 - errors_10;

            let chr = &line[0];
            let start = &line[1].parse::<i32>().unwrap();
            let end = &line[2].parse::<i32>().unwrap();
            println!("{}\t{}\t{}\t{}\t{}\t{}\n", chr, start, end, t, m, cell_type);
        }
    }
    Ok(())
}
#[pymodule]
fn _nanomix(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<MMSE>()?;
    m.add_function(wrap_pyfunction!(generate_methylome, m)?)?;
    Ok(())
}


#[cfg(test)]
mod tests {
    use std::io::prelude::*;
    use std::fs::File;
    use crate::*;
    use ndarray::arr2;
    #[test]
    fn init_bmm() {
        let mut buffer = File::create("methylome.test.tsv").unwrap();
        buffer.write_all(b"chromosome\tstart\tend\tmodified_calls\ttotal_calls\n").unwrap();
        buffer.write_all(b"chr1\t10\t200\t5\t10\n").unwrap();
        buffer.write_all(b"chr1\t205\t500\t1\t3\n").unwrap();
        let mut buffer2 = File::create("atlas.test.tsv").unwrap();
        buffer2.write_all(b"chromosome\tstart\tend\tadipocytes\tmonocytes\n").unwrap();
        buffer2.write_all(b"chr1\t5\t510\t0.9\t0.3\n").unwrap();
        buffer2.write_all(b"chr2\t5\t510\t0.5\t0.4\n").unwrap();
        let model = MMSE::new("methylome.test.tsv", "atlas.test.tsv", 0.05, 0.95).unwrap();
        assert_eq!(model.A, arr2(&[[0.9, 0.3], [0.9, 0.3]]));
        assert_eq!(model.m, vec![5, 1]);
        assert_eq!(model.t, vec![10, 3]);
        assert_eq!(model.sigma, vec![0.5, 0.5]);
    }
    //#[test]
    fn init_bmm_with_generated_data() {
        generate_methylome("adipocyte_simulated_methylome.tsv", "loyfer250Atlas.tsv", 0.05, 0.95, 0);
        let model = MMSE::new("adipocyte_simulated_methylome.tsv", "loyfer250Atlas.tsv", 0.05, 0.95).unwrap();
    }
    #[test]
    fn gamma_tilde() {
        let model = MMSE::new("methylome.test.tsv", "atlas.test.tsv", 0.0, 1.0).unwrap();
        assert!(model.gamma_tilde(0, 0) - 0.5*0.00149 < 0.00001, "Left = {}, Right = {}", model.gamma_tilde(0, 0), 0.5*0.00149);
        assert!(model.gamma_tilde(1, 0) - 0.5*0.10292 < 0.00001, "Left = {}, Right = {}", model.gamma_tilde(0, 0), 0.5*0.10292);
    }
    #[test]
    fn gamma() {
        let model = MMSE::new("methylome.test.tsv", "atlas.test.tsv", 0.0, 1.0).unwrap();
        assert_eq!(model.gamma(0, 0) + model.gamma(1, 0), 1.0)
    }
    #[test]
    fn optimize() {
        let mut model = MMSE::new("methylome.test.tsv", "atlas.test.tsv", 0.0, 1.0).unwrap();
        model.optimize(5.0, 50);
        assert_eq!(model.sigma.iter().sum::<f64>(), 1.0)
    }
}
