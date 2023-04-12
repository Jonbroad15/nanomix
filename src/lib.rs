use pyo3::prelude::*;
use rand::distributions::{Weighted, WeightedChoice, Distribution, Poisson};
use rand::prelude::*;
use csv;
use std::io::{BufReader, BufRead};
use std::fs;
use intervaltree::*;
use std::collections::HashMap;
use core::ops::Range;
use rulinalg::vector::Vector;
use logaddexp::{LogSumExp};

#[derive(Debug)]
pub struct Read {
    name: String,
    chr: String,
    start: usize,
    end: usize,
    atlas_idx: usize,
    m: usize,
    t: usize,
}
#[derive(Clone, Debug)]
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
            let values: Vec<f64> = record.iter().skip(3).map(|x| x.parse().unwrap_or(0.0)).collect();

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
        // print the atlas shape
        eprintln!("Atlas shape: {} x {}", atlas_data.len(), atlas_data[0].len());
        Atlas {
            region_map: interval_trees,
            values: atlas_data,
            cell_types: cell_types,
        }
    }
    pub fn query(&self, chr: &String, start: usize, end: usize) -> Result<usize, &str> {
        //Query the region_map to get the methylation propensity
        // check if the read.chr is in the region_map
        if let Some(tree) = self.region_map.get(chr) {
            let mut query = tree.query(start..end);
            // if there is a match, get the index of the interval
            // if there is no match, return an error
            // if there are multiple matches, return an error
            if let Some(element) = query.next() {
                if query.next().is_some() {
                    return Err("Multiple matches found");
                }
                Ok(element.value)
            } else {
                Err("No match found")
            }
        } else {
            return Err("Chromosome not found in region map")
        }
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
    p01: Vec<f64>,
    p11: Vec<f64>,
    read_name_map: HashMap<String, Vec<usize>>,
    #[pyo3(get)]
    pub assignments: Vec<usize>,
    concentration: f64,
    gamma_sum: HashMap<String, f64>
}
impl MMSE {
    fn get_p(&self, xi: f64, k: usize) -> f64 {
        self.p11[k]*xi + self.p01[k]*(1.0 - xi)
    }
    fn gamma_tilde(&self, k: usize, p: f64, read: &Read) -> Result<f64, &str> {
        //let bin = Binomial::new(read.t, p);
        //compute the binomial coefficient for read.t choose read.m
        // clip p to avoid numerical issues
        let p = p.min(1.0 - 1e-10);
        let p = p.max(1e-10);
        let logpmf: f64 = read.m as f64 * p.ln() + (read.t - read.m) as f64 * (1.0 - p).ln();
        let gamma_tilde = self.sigma[k] + logpmf; //+ ln_binomial(read.t as u64, read.m as u64);
        if gamma_tilde.is_nan() {
            eprintln!("k = {}, self.sigma[k] = {}, p = {}, read = {:?}\n", k, self.sigma[k], p, read);
            eprintln!("self.sigma = {:?}\n", self.sigma);
            eprintln!("logpmf = {}\n", logpmf);
            return Err("gamma_tilde is NaN")
        } else {
            return Ok(gamma_tilde)
        }
    }
    fn gamma_tilde_joint(&self, k: usize, read: &String) -> Result<f64, &str> {
        // loop through all the Reads assoicated with the read name
        // compute the joint probability of the reads
        let mut gamma_tilde_joint = 0.0;
        for read in self.read_name_map.get(read).unwrap() {
            let read = &self.methylome[*read];
            let p = self.get_p(self.atlas.values[read.atlas_idx][k], k);
            gamma_tilde_joint += self.gamma_tilde(k, p, read)?;
        }
        Ok(gamma_tilde_joint)
    }
    pub fn gamma(&self, k: usize, read: &String) -> Result<f64, &str> {
        let gamma: f64 = self.gamma_tilde_joint(k, read).unwrap() - self.gamma_sum.get(read).unwrap();

        if gamma.is_nan()  {
            eprintln!("self.gamma_tilde_joint(k, self.get_p(atlas_row[k]), read, &self.sigma).unwrap() = {}\n",
            self.gamma_tilde_joint(k, read).unwrap());
            eprintln!("k = {}, sigma[k] = {}, read = {}\n", k, self.sigma[k], read);
            eprintln!("sigma = {:?}\n", self.sigma);
            eprintln!("sum = {}\n", self.gamma_sum.get(read).unwrap());
            eprintln!("gamma = {}\n", gamma);
            return Err("gamma is NaN")
        } else {
            return Ok(gamma)
        }
    }
    fn print_gamma(&self, i: usize) -> (){
        let read = self.read_name_map.keys().nth(i).unwrap();
        (0..self.K).for_each(|k| eprint!("{}:{:.1}\t",k,self.gamma(k, read).unwrap()));
        eprintln!("gamma_sum: {}",
                 (0..self.K).map(|k| self.gamma(k, read).unwrap()).ln_sum_exp());
    }
    fn print_sigma(&self) -> (){
        let sum = (0..self.K).map(|k| self.sigma[k].exp()).sum::<f64>();
        (0..self.K).for_each(|k| eprint!("{}:{:.2}\t", k , self.sigma[k].exp() / sum));
        eprintln!("sigma_sum: {}", sum);
    }
    fn print_gamma_tilde(&self, i: usize) -> (){
        // Print the gamma_tilde_joint values for a given read
        let read = self.read_name_map.keys().nth(i).unwrap();
        (0..self.K).for_each(|k| eprint!("{}:{:.1}\t",k,self.gamma_tilde_joint(k, read,).unwrap()));
        eprintln!("gamma_tidle_sum: {}",
                 (0..self.K).map(|k| self.gamma_tilde_joint(k, read).unwrap()).ln_sum_exp());
    }
    fn e_step(&mut self) -> f64{
        // compute the log likelihood of the data
        let mut acc = 0.0;
        let mut gamma_sum = 0.0;
        for read in self.read_name_map.keys() {
            let mut gammas = vec![];
            for k in 0..self.K {
                gammas.push(self.gamma_tilde_joint(k, read).unwrap());
            }
            gamma_sum = gammas.into_iter().ln_sum_exp();
            // update the gamma_sum map
            self.gamma_sum.insert(read.to_string(), gamma_sum);
            acc += gamma_sum;
        }
        acc
    }
    fn remove_uninformative_reads(&mut self, t: f64) -> () {
        // remove reads such that their gamma values are equal to sigma for every cell type k
        // this is a heuristic to remove reads that are not informative
        let mut to_remove = vec![];
        for (i, read) in self.read_name_map.keys().enumerate() {
            let mut remove = true;
            for k in 0..self.K {
                if (self.gamma(k, read).unwrap().exp()  - self.sigma[k].exp()).abs() > t {
                    remove = false;
                    break;
                }
            }
            if remove {
                to_remove.push(read.clone());
            }
        }
        for i in to_remove {
            self.read_name_map.remove(&i);
        }
        // update N
        self.N = self.read_name_map.len();
    }
    fn m_step(&mut self){
        let mut sigma = vec![0.0; self.K];
        for k in 0..self.K {
            let mut gammas = vec![];
            for read in self.read_name_map.keys() {
                gammas.push(self.gamma(k, read).unwrap());
            }
             //if all the gammas are -inf, then set sigma to -inf
            //if gammas.iter().all(|x| x.is_infinite() && x.is_sign_negative()) {
                //sigma[k] = 0.0;
            //} else {
            sigma[k] = gammas.clone().into_iter().map(|x| x.exp() ).sum::<f64>() / (self.N as f64);
            //}
            if sigma[k].is_nan() {
                eprintln!("gammas = {:?}\n", gammas);
                eprintln!("k = {}, self.sigma[k] = {}\n", k, self.sigma[k]);
                eprintln!("self.sigma = {:?}\n", self.sigma);
                eprintln!("self.N.ln() = {}\n", (self.N as f64).ln());
                eprintln!("gammas.into_iter().ln_sum_exp() = {}\n", 
                          gammas.clone().into_iter().ln_sum_exp());
                panic!("sigma is NaN")
            }
        }
        // adjust so sigma adds up to 1
        //let sigma_sum = sigma.iter().sum::<f64>();
        //sigma.iter_mut().for_each(|x| *x = *x / sigma_sum);
        // convert to log space
        self.sigma = sigma.into_iter().map(|x| x.ln()).collect();
    }
    fn noise_trim(&mut self, t: f64) {
        // take sigma out of ln space
        let mut sigma = self.sigma.iter().map(|x| x.exp()).collect::<Vec<f64>>();
        let mut sum: f64 = sigma.iter().sum();
        sigma = sigma.iter().map(|x| x/sum).collect::<Vec<f64>>();
        for k in 0..self.K {
            if sigma[k] < t { sigma[k] = 0.0 }
        }
        sum = sigma.iter().sum();
        // if sum is 0, set all sigmas to 1/K
        if sum == 0.0 {
            eprintln!("Noise trim set all sigmas to 1/K, since no proportions were above threshold");
            self.sigma = (0..self.K).map(|_| (1.0/self.K as f64).ln()).collect();
        } else {
            self.sigma = sigma.iter().map(|x| (x/sum).ln()).collect();
        }
    }
}


#[pymethods]
impl MMSE {
    #[new]
    pub fn new(methylome: &str, atlas: &str, init_sigma: Vec<f64>, p01: Vec<f64>, p11: Vec<f64>, concentration: f64) -> PyResult<Self> {
        let mut atlas = Atlas::new(atlas);
        // Read methylome into memory
        // Store the list of intervals for the methylome instead of interval tree
        // Make a struct for reads
        let mut methylome_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(methylome)
            .expect("Could not open methylome file");

        // get header from methylome reader
        let header = methylome_reader.headers().unwrap();
        // Map header to index
        let header_map = header.iter().enumerate().map(|(i, x)| (x.to_string(), i)).collect::<HashMap<String, usize>>();
        // get the index of the columns we need
        let read_name_idx = header_map.get("read_name").unwrap();
        let chr_idx = header_map.get("chromosome").unwrap();
        let start_idx = header_map.get("start_position").unwrap();
        let end_idx = header_map.get("end_position").unwrap();
        let m_idx = header_map.get("modified_calls").unwrap();
        let t_idx = header_map.get("total_calls").unwrap();
        // this stores the methylome data for each interval
        let mut methylome_data = Vec::new();
        // enumerate methylome_reader.records()
        for (i, r) in methylome_reader.records().enumerate() {
            let record = r.expect("Could not parse bed record");
            let chr = record[*chr_idx].to_string();
            let start = record[*start_idx].parse::<usize>().expect(&format!("Could not parse start position at index: {}", i));
            let end = record[*end_idx].parse::<usize>().expect(&format!("Could not parse end position at index: {}", i));
            // check if read is in the atlas region
            if let Ok(idx) = atlas.query(&chr, start, end) {
                    let read = Read {
                        chr: chr,
                        start: start,
                        end: end,
                        name: record[*read_name_idx].to_string(),
                        m: record[*m_idx].parse::<usize>().expect(&format!("Could not parse modified calls at index: {}", i)),
                        t: record[*t_idx].parse::<usize>().expect(&format!("Could not parse total calls at index: {}", i)),
                        atlas_idx: idx,
                    };  
                    methylome_data.push(read);
                }
        };
        let K = atlas.values[0].len();




        // make a hashmap of read names to index in methylome_data
        // keys are read names and values are vectors of indices
        let mut read_name_map = HashMap::new();
        for (i, read) in methylome_data.iter().enumerate() {
            // if the read name is not in the hashmap, add it
            // if the read name is in the hashmap, make a new read name with a number at the end
            // and add that to the hashmap
            let mut read_name = read.name.clone();
            let mut j = 0;
            while read_name_map.contains_key(&read_name) {
                read_name = format!("{}_{}", read.name, j);
                j += 1;
            }
            read_name_map.entry(read_name.clone()).or_insert(vec![]).push(i);
        }
        let mut gamma_sum = HashMap::<String, f64>::new();
        for read in methylome_data.iter() {
             gamma_sum.insert(read.name.clone(), 0.0);
        }
        // calculate N to be the number of keys in the read_name_map
        let N = read_name_map.keys().len();
        eprintln!("total methylated calls: {}", methylome_data.iter().map(|x| x.m).sum::<usize>());
        eprintln!("total unmethylated calls: {}", methylome_data.iter().map(|x| x.t - x.m).sum::<usize>());
        eprintln!("total reads: {}", N);

        Ok(MMSE{
            sigma: init_sigma.iter().map(|x| x.ln()).collect(),
            N: N,
            K: K,
            atlas: atlas,
            methylome: methylome_data,
            p01: p01,
            p11: p11,
            read_name_map: read_name_map,
            assignments: vec![0 as usize; N],
            concentration: concentration,
            gamma_sum: gamma_sum})
    }
    pub fn optimize(&mut self, stop_thresh: f64, max_iter: u32, min_proportion: f64)-> PyResult<()>{
        let mut i = 0;
        let mut ll_prev = self.e_step();
        let mut ll = ll_prev;
        
        eprintln!("Running EM Inference");
        loop {
            i += 1;
            self.m_step();
            self.print_sigma();
            //if i == 10 {
                //eprintln!("Removing uninformative regions, N = {}", self.N);
                //self.remove_uninformative_reads(0.01);
                //self.trim_noise(min_proportion);
                //eprintln!("After removal N = {}", self.N);
            //}
            ll = self.e_step();
            let condition = 100.0*((ll - ll_prev)/ll_prev.abs());
            eprintln!("Iteration: {}\tLog-likelihood: {:.5}\tPercent_change: {:.8}",
                      i, ll, condition);
            if ( condition.abs() < stop_thresh && i > 10 ) || i >= max_iter-1 {
                break;
            } else if  condition < 0.0 {
                eprintln!("Warning: Log-likelihood decreased");
                let cell_type_proportions = self.cell_type_proportions();
                //eprintln!("Cell type proportions: {:?}", cell_type_proportions);
                //print gamma for the first 10 reads
                self.print_gamma(0);
                self.print_gamma_tilde(0);
                ll_prev = ll;
            } else {
                ll_prev = ll;
            }
        }
        eprintln!("Sigma: {:?}", self.sigma);
        self.noise_trim(min_proportion);
        eprintln!("Sigma: {:?}", self.sigma);
        self.e_step();
        self.m_step();
        eprintln!("Sigma: {:?}", self.sigma);
        ll = self.e_step();
        i += 1;
        eprintln!("Model optimized in {} steps. Log-likelihood = {}", i, ll);
        eprintln!("log-likelihood:\t{:.2}", ll);
        eprintln!(" N = {}", self.N);
        Ok(())
    }
    pub fn evaluate(&mut self, stop_thresh: f64, max_iter: u32, min_proportion: f64,
                    sigma: Vec<f64>, target_assignments: Vec<usize>) -> PyResult<()>{
        let mut i = 0;
        let mut ll_prev = self.e_step();
        let mut ll = self.e_step();
        loop {
            let loss = deconvolution_loss(&sigma, &self.sigma);
            let acc = self.accuracy(target_assignments.clone(), 0.5);
            eprintln!("Iter: {}\tLoss: {:.3}\tAccuracy: {:.3}\tLog-Likelihood: {:.3}\t True Log-likelihood: {:.3}",
                      i, loss, acc, ll, self.e_step());
            i += 1;

            self.m_step();
            ll = self.e_step();
            if (ll_prev - ll).abs()/ll < stop_thresh || i == max_iter {
                break;
            } else {
                ll_prev = ll;
            }
        }
        self.noise_trim(min_proportion);
        self.e_step();
        self.m_step();
        ll = self.e_step();
        let loss = deconvolution_loss(&sigma, &self.sigma);
        let acc = self.accuracy(target_assignments.clone(), 0.5);
        eprintln!("Loss: {:.3}\t Accuracy: {:.3}\tLog-Likelihood: {:.3}\t True Log-likelihood: {:.3}", loss, acc, ll, self.e_step());
        Ok(())
    }
    pub fn reset(&mut self){
        self.sigma = (0..self.K).map(|_k| 1.0/self.K as f64).collect()
    }
    fn assign_fragments(&mut self) -> () {
        for (i, read) in self.read_name_map.keys().enumerate() {
            let gamma_vec: Vec<f64> = (0..self.K).map(|k| self.gamma(k, read).unwrap()).collect();
            self.assignments[i] = Vector::new(gamma_vec).argmax().0;
        }
        ()
    }
    pub fn assign_fragments_t(&mut self, threshold: f64) -> Vec<usize> {
        // assign reads to the most likely cell type with confidence above threshold
        // If the read is unnassigned it gets value K
        let mut assignments = vec![0 as usize; self.N];
        self.e_step();
        for (i, read) in self.read_name_map.keys().enumerate() {
            let gamma_vec: Vec<f64> = (0..self.K).map(|k| self.gamma(k, read).unwrap()).collect();
            if gamma_vec.iter().cloned().fold(0./0., f64::max) > threshold.ln() {
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
    pub fn log_likelihood(&mut self) -> f64 {
        self.e_step()
    }
    pub fn cell_type_proportions(&mut self) -> HashMap<String, f64> {
        // Return a hashmap of keys = cell type names and values = respective sigma proportions
        let mut cell_type_proportions : HashMap<String, f64> = HashMap::new();
        // Convert sigma from log space and normalize to sum to 1
        let sigma = self.sigma.iter().map(|x| x.exp()).collect::<Vec<f64>>();
        let sigma_sum = sigma.iter().sum::<f64>();
        let sigma = sigma.iter().map(|x| x/sigma_sum).collect::<Vec<f64>>();
        for (k, cell_type) in self.atlas.cell_types.iter().enumerate() {
            cell_type_proportions.insert(cell_type.clone(), sigma[k]);
        }
        cell_type_proportions
    }
    //pub fn class_probabilities(&mut self, target_assignments: Vec<usize>) -> Vec<Vec<f64>> {
        //self.assign_fragments();
        //let mut probs: Vec<Vec<f64>> = self.assignments.iter().zip(target_assignments.iter())
            //.map(|(value, target)| (*value == *target) as u32 as f64 )
            //.map(|correct| vec![correct] )
            //.collect();
        //(0..self.N).map(|i| (self.atlas.query(&self.methylome[i]).unwrap().clone(), &self.methylome[i]) )
            //.map(|(atlas_row, read)| (0..self.K).map(|k| self.gamma(k, atlas_row, read).unwrap()).collect())
            //.enumerate()
            //.for_each(|(i, mut row)| probs[i].append(&mut row));
        //(0..self.N)
            //.map(|i| self.atlas.query(&self.methylome[i]).unwrap().clone() ).enumerate()
            //.for_each(|(i, mut row)| probs[i].append(&mut row) );
        //probs
    //}
}
pub fn deconvolution_loss(sigma: &Vec<f64>, sigma_hat: &Vec<f64>) -> f64 {
    sigma.iter().zip(sigma_hat.iter()).map(|(y, yhat)| (y-yhat).powf(2.0)).sum::<f64>().powf(0.5)
}
#[pyfunction]
pub fn generate_methylome(atlas: &str, sigma: Vec<f64>, coverage: f64, region_size: u64, p01: f64, p11: f64) -> PyResult<()> {
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
    if cell_types.len() != sigma.len() { panic!("Length of Sigma ({}) != number of cell types ({})", sigma.len(), cell_types.len()) }

    let poi = Poisson::new(coverage);
    let mut rng = StdRng::from_entropy();
    let mut choices: Vec<Weighted<usize>> = sigma.iter().enumerate()
        .map(|(i, x)| Weighted { weight: (x*1000.0) as u32, item: i })
        .collect();
    let wc = WeightedChoice::new(&mut choices);

    println!("read_name\tchromosome\tstart_position\tend_position\ttotal_calls\tmodified_calls\tcell_type\n");

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
            // generate a random read name
            //let read_name = rand::thread_rng().sample_iter(&Alphanumeric).take(10).collect::<String>();
            // generate a random read name unique to chr start end
            let read_name = format!("{}_{}_{}_{}", chr, start, end, coverage);
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", read_name, chr, start, end, t, m, cell_type);
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
    #[test]
    fn test_cell_type_proportions() {
        let mut buffer = File::create("methylome.test.tsv").unwrap();
        buffer.write_all(b"chromosome\tstart\tend\tmodified_calls\ttotal_calls\n").unwrap();
        buffer.write_all(b"chr1\t10\t200\t5\t10\n").unwrap();
        buffer.write_all(b"chr1\t205\t500\t1\t3\n").unwrap();
        let mut buffer2 = File::create("atlas.test.tsv").unwrap();
        buffer2.write_all(b"chromosome\tstart\tend\tadipocytes\tmonocytes\n").unwrap();
        buffer2.write_all(b"chr1\t5\t510\t0.9\t0.3\n").unwrap();
        buffer2.write_all(b"chr2\t5\t510\t0.5\t0.4\n").unwrap();
        let sigma = vec![0.7, 0.3];
        let mut model = MMSE::new("methylome.test.tsv", "atlas.test.tsv", sigma, 0.05, 0.95).unwrap();
        let cell_type_proportions = model.cell_type_proportions();
        assert_eq!(cell_type_proportions, HashMap::from_iter(vec![("adipocytes".to_string(), 0.7), ("monocytes".to_string(), 0.3)]));
    }
}
