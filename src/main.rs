use clap::{command, Parser};
use csv::{Reader, Writer};
use num::{Float, NumCast};
use std::fs::File;
use std::iter::Iterator;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    author = "Luca Visentin",
    about = "Calculate cohen's d of expression values."
)]
struct Args {
    /// Path to the expression matrix with the 'case' expression matrix
    case_expression_matrix: PathBuf,
    /// Path to the expression matrix with the 'control' expression matrix
    control_expression_matrix: PathBuf,
    /// Path and filename of the output file
    output_path: PathBuf,
    /// Delimiter of the input files
    #[arg(short, long, default_value_t = ',')]
    delimiter: char,
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

    let mut writer =
        Writer::from_path(&args.output_path).expect("Could not open output file for writing.");

    process_csvs(&mut case_samples, &mut control_samples, &mut writer);
}

fn process_csvs<R, T, V>(
    case_samples: &mut Reader<R>,
    control_samples: &mut Reader<T>,
    writer: &mut Writer<V>,
) where
    V: std::io::Write,
    R: std::io::Read + std::io::Seek,
    T: std::io::Read + std::io::Seek,
{
    let control_start_position = control_samples.position().clone();
    let case_start_position = case_samples.position().clone();
    // Sort the control and case samples to have the same row names
    // I assume that the first row is made up of column names
    let control_row_names: Vec<String> = control_samples
        .records()
        .map(|x| x.unwrap().get(0).unwrap().to_owned())
        .collect();
    let case_row_names: Vec<String> = case_samples
        .records()
        .map(|x| x.unwrap().get(0).unwrap().to_owned())
        .collect();

    case_samples.seek(case_start_position).unwrap();
    control_samples.seek(control_start_position).unwrap();

    // I clone since I have to re-borrow this later...
    let row_names_match = control_row_names
        .clone()
        .into_iter()
        .zip(case_row_names)
        .all(|(x, y)| x == y);

    if !row_names_match {
        println!("ERROR: Row names between case and control files do not match up");
        return;
    };

    println!("Computing cohen's d...");
    let result: Vec<f64> = case_samples
        .records()
        .zip(control_samples.records())
        .skip(1)
        .map(|(case, control)| {
            let case_values: Vec<f64> = case
                .expect("Couldn't read case record")
                .into_iter()
                .skip(1)
                .map(|x| {
                    x.parse()
                        .unwrap_or_else(|_| panic!("non-float value in case record: {x}"))
                })
                .collect();

            let control_values: Vec<f64> = control
                .expect("Couldn't read case record")
                .into_iter()
                .skip(1)
                .map(|x| {
                    x.parse()
                        .unwrap_or_else(|_| panic!("non-float value in case record: {x}"))
                })
                .collect();

            cohen(case_values, control_values)
        })
        .collect();

    writer.write_record(vec!["row_names", "cohen_d"]).unwrap();
    for (row_name, value) in control_row_names.into_iter().zip(result.into_iter()) {
        writer
            .write_record(vec![row_name, format!("{}", value)])
            .unwrap();
    }
    writer.flush().unwrap();

    println!("Done!");
}

/// Calculate the mean of the values in the input vector
fn mean<F>(data: Vec<F>) -> Option<F>
where
    F: Float + std::iter::Sum,
{
    let count = data.len();
    let sum = data.into_iter().sum::<F>();

    match count {
        positive if positive > 0 => Some(sum / NumCast::from(count).unwrap()),
        _ => None,
    }
}

/// Calculate the variance of the values in the input vector
fn var<F>(data: &Vec<F>) -> F
where
    F: Float + std::iter::Sum,
{
    let data_mean = mean(data.clone()).unwrap();
    let count = &data.len();

    let variance = data
        .iter()
        .map(|value| {
            let diff = (*value as F) - data_mean;

            diff * diff
        })
        .sum::<F>()
        / NumCast::from(*count - 1).unwrap();

    variance
}

/// Calculate cohen's d statistic from a case and control numeric vectors.
fn cohen<T>(case: Vec<T>, control: Vec<T>) -> T
where
    T: Float + std::iter::Sum,
{
    let n_case: T = NumCast::from(case.len()).unwrap();
    let n_control: T = NumCast::from(control.len()).unwrap();

    let pooled_stdev = (((n_case - NumCast::from(1).unwrap()) * var(&case)
        + (n_control - NumCast::from(1).unwrap()) * var(&control))
        / (n_case + n_control - NumCast::from(2).unwrap()))
    .powf(NumCast::from(0.5).unwrap());

    if pooled_stdev == NumCast::from(0).unwrap() {
        return NumCast::from(0).unwrap();
    }

    (mean(case).unwrap() - mean(control).unwrap()) / pooled_stdev
}

#[cfg(test)]
mod tests {
    use crate::*;
    use csv::{ReaderBuilder, WriterBuilder};
    use num::{abs, Signed};
    use std::io::Cursor;
    use std::ops::Sub;

    fn kinda_equal<T, G>(a: T, b: T, tolerance: G) -> bool
    where
        T: Sub + PartialEq,
        G: PartialEq,
        <T as Sub>::Output: PartialOrd<G>,
        <T as Sub>::Output: Signed,
    {
        abs(a - b) < tolerance
    }

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

    #[test]
    fn test_cohen() {
        assert!(kinda_equal(
            cohen(vec![2.2, 1.3, 3.1], vec![12.6, 11.1, 12.3, 13.2, 14.2]),
            -9.809684, // Calculated this by hand
            0.0001
        ))
    }

    #[test]
    fn itegration() {
        let cases = "\
row_names,sample1,sample2,sample3
gene_1,2.2,1.3,3.1
gene_2,1.3,2.2,3.1
gene_3,3.1,2.2,1.3
";
        let controls = "\
row_names,sample4,sample5,sample6
gene_1,12.6,11.1,12.3
gene_2,11.1,12.3,12.6
gene_3,12.3,12.6,11.1
";
        let mut case_samples = ReaderBuilder::new().from_reader(Cursor::new(cases.as_bytes()));
        let mut control_samples =
            ReaderBuilder::new().from_reader(Cursor::new(controls.as_bytes()));
        let mut writer = WriterBuilder::new().from_writer(vec![]);

        process_csvs(&mut case_samples, &mut control_samples, &mut writer);

        let output_data = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output_data,
            "\
row_names,cohen_d
gene_1,-11.549410759380276
gene_2,-11.549410759380276
gene_3,-11.549410759380276
"
        )
    }

    #[test]
    fn itegration_real() {
        let cases = "\
row_names,sample1,sample2,sample3
gene_1,-0.535058383960151,1.36337028207967,1.94555121778008
gene_2,-0.00508937611274737,2.05525510002174,-0.605214688134933
gene_3,-0.0224148926648615,1.00937541977038,-0.675006944468402
";
        let controls = "\
row_names,sample4,sample5,sample6,sample7
gene_1,2.6842520940552,5.2521840950163,3.75921244478561,4.83016238235602
gene_2,-2.29756231844181,-1.15137785595405,-0.486929672999351,-2.816582037462001
gene_3,9.6425133676228,11.0997575073032,9.17697194351323,8.82739260182
";
        let mut case_samples = ReaderBuilder::new().from_reader(Cursor::new(cases.as_bytes()));
        let mut control_samples =
            ReaderBuilder::new().from_reader(Cursor::new(controls.as_bytes()));
        let mut writer = WriterBuilder::new().from_writer(vec![]);

        process_csvs(&mut case_samples, &mut control_samples, &mut writer);

        let output_data = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output_data,
            "\
row_names,cohen_d
gene_1,-2.6462737190827195
gene_2,1.7993829062243292
gene_3,-10.169673185162345
"
        )
    }
}
