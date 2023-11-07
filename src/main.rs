use clap::{command, Parser};
use csv::Reader;
use std::fs::File;
use std::path::PathBuf;
use std::iter::Iterator;
use num::{Float, NumCast};

#[derive(Parser, Debug)]
#[command(
    author = "Luca Visentin",
    about = "Calculate cohen's D of expression values."
)]
struct Args {
    /// Path to the expression matrix with the 'case' expression matrix
    case_expression_matrix: PathBuf,
    /// Path to the expression matrix with the 'control' expression matrix
    control_expression_matrix: PathBuf,
    /// Path and filename of the output file
    output_path: PathBuf,
    /// Delimiter of the input files
    #[arg(short, long, default_value_t = '\t')]
    delimiter: char
}

fn main() {
    // Parse the command-line arguments.
    let args: Args = Args::parse();
    
    println!("Building case reader...");
    let mut case_samples: Reader<File> = csv::ReaderBuilder::new()
        .delimiter(args.delimiter as u8)
        .from_path(&args.case_expression_matrix)
        .expect("Failed to read input case expression matrix");


    println!("Building control reader...");
    let mut control_samples: Reader<File> = csv::ReaderBuilder::new()
        .delimiter(args.delimiter as u8)
        .from_path(&args.control_expression_matrix)
        .expect("Failed to read input case expression matrix");

    // Sort the control and case samples to have the same row names
    // I assume that the first row is made up of column names
    let control_row_names: Vec<String> = control_samples.records()
        .map(|x| x.unwrap().get(0).unwrap().to_owned()).collect();
    let case_row_names: Vec<String> = case_samples.records()
        .map(|x| x.unwrap().get(0).unwrap().to_owned()).collect();

    let row_names_match = control_row_names.into_iter().zip(case_row_names.into_iter()
    ).all(|(x, y)| x == y);

    if !row_names_match {
        println!("ERROR: Row names between case and control files do not match up.");
        return ();
    };

    println!("Computing cohen's D...");
    let _result: Vec<f64> = case_samples.records().zip(control_samples.records()).map(
        |(case, control)| {
            let case_values: Vec<f64> = case
                .expect("Couldn't read case record")
                .into_iter()
                .skip(1)
                .map(|x| x.parse().expect(&format!("non-float value in case record: {x}")))
                .collect();
            
            let control_values: Vec<f64> = control
                .expect("Couldn't read case record")
                .into_iter()
                .skip(1)
                .map(|x| x.parse().expect(&format!("non-float value in case record: {x}")))
                .collect();


            cohen(case_values, control_values)
        }
    ).collect();

    ()
}

/// TODO: The mean, var and cohen functions would be more useful if we made them
/// generic. But I don't know how to do that (very well), so I leave them be.

/// Calculate the mean of the values in the input vector
fn mean<F>(data: Vec<F>) -> Option<F> where F: Float + std::iter::Sum {
    let count = data.len();
    let sum = data.into_iter().sum::<F>();

    match count {
        positive if positive > 0 => Some(sum / NumCast::from(count).unwrap()),
        _ => None,
    }
}

/// Calculate the variance of the values in the input vector
fn var<F>(data: &Vec<F>) -> F where F: Float + std::iter::Sum {
    let data_mean = mean(data.clone()).unwrap();
    let count = &data.len();

    let variance = data
        .iter()
        .map(|value| {
            let diff = data_mean - (*value as F);

            diff * diff
        })
        .sum::<F>()
        / NumCast::from(*count - 1).unwrap();

    variance
}

/// Calculate cohen's D statistic from a case and control numeric vectors.
fn cohen<T>(case: Vec<T>, control: Vec<T>) -> T where T: Float + std::iter::Sum {
    let n_case: T = NumCast::from(case.len()).unwrap();
    let n_control: T = NumCast::from(case.len()).unwrap();

    let pooled_var = (
        ((n_case - NumCast::from(1).unwrap()) * var(&case) + (n_control - NumCast::from(1).unwrap()) * var(&control)) / (n_case + n_control - NumCast::from(2).unwrap())
    ).powf(NumCast::from(0.5).unwrap());

    if pooled_var == NumCast::from(0).unwrap() {
        return NumCast::from(0).unwrap();
    }

    (mean(case).unwrap() - mean(control).unwrap()) / pooled_var
    //(mean(case).unwrap() - mean(control).unwrap()) as f32
}


#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn mean_of_values() {
        assert_eq!(mean(vec![1., 2., 3.]).unwrap(), 2.);
        assert_eq!(mean(vec![10., 10., 10.]).unwrap(), 10.);
        assert_eq!(mean(vec![0., 12., 0., 23.]).unwrap(), 8.75);
        assert_eq!(mean(vec![1., 2., 1.]).unwrap(), 1.3333333333333333);
    }

    #[test]
    fn variance_of_values() {
        assert_eq!(var(&vec![1., 2., 3.]), 1.);
        assert_eq!(var(&vec![10., 10., 10.]), 0.);
        assert_eq!(var(&vec![0., 12., 0., 23.]), 122.25);
        assert_eq!(var(&vec![1., 2., 1.]), 0.3333333333333333);
    }
}