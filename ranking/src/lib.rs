use ndarray::{Axis, Array2, ArrayView2, ArrayView1};
use numpy::{IntoPyArray, PyReadonlyArray2, PyArray2};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};
use std::collections::HashMap;


#[pymodule]
fn aln_ranking(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    fn rank_sequences(seq: ArrayView2<'_, u32>) -> Array2<u32> {

        let x: Vec<_> = seq.axis_iter(Axis(1))
            .map(|col| _count_occurrences(&col))
            .collect();
        let y: Vec<_> = x.iter()
            .map(|occ| _sort_occurrences(&occ))
            .collect();
        let z: Vec<_> = y.iter()
            .map(|occ| _convert_to_ranks(&occ))
            .collect();

        let mut seq_ranked: Array2<u32> = Array2::from_elem(seq.raw_dim(), 0);
        for (i, col) in seq.axis_iter(Axis(1)).enumerate() {
            let current_ranks = &z[i];
            for (j, elem) in col.iter().enumerate() {
                seq_ranked[[j, i]] = current_ranks[elem];
            }
        }
        seq_ranked
        }

        #[pyfn(m, "rank_sequences")]
        fn rank_sequences_py<'py>(
            py: Python<'py>, seq: PyReadonlyArray2<'_, u32>
        ) -> &'py PyArray2<u32> {
            rank_sequences(seq.as_array()).into_pyarray(py)
        }

        Ok(())
    }


fn _count_occurrences(seq: &ArrayView1<u32>) -> HashMap<u32, u32> {
    let mut char_counts: HashMap<u32, u32> = HashMap::new();

    for &c in seq {
        if c != 0 {
            *char_counts.entry(c).or_insert(0) += 1;
        }
    }
    char_counts
}

fn _sort_occurrences(occur_map: &HashMap<u32, u32>) -> Vec<(&u32, &u32)> {
    let mut count_vec: Vec<_> = occur_map.iter().collect();
    count_vec.sort_by(|a, b| {
        b.1.cmp(a.1).then(a.0.cmp(b.0))
    });
    count_vec
}

fn _convert_to_ranks(sorted_occurrs: &Vec<(&u32, &u32)>) -> HashMap<u32, u32> {
    let mut ranks: HashMap<u32, u32> = HashMap::new();
    let occurs_len = sorted_occurrs.len();
    for (i, pair) in sorted_occurrs.iter().enumerate()  {
        ranks.insert(*pair.0, (occurs_len - i) as u32);
    };
    ranks.insert(0, 0);
    ranks
}

mod tests {
    use super::*;

    #[test]
    fn test_count_occurences() {
        let input = ArrayView1::from(&[0, 2, 0, 1, 0, 2, 2]);

        let obs = _count_occurrences(&input);
        let mut exp = HashMap::new();
        exp.insert(1, 1);
        exp.insert(2, 3);

        assert_eq!(obs, exp)
    }

    #[test]
    fn test_sort_occurences() {
        let mut input = HashMap::new();
        input.insert(1, 2);
        input.insert(2, 3);
        input.insert(3, 4);
        input.insert(4, 4);
        input.insert(6, 2);

        let obs = _sort_occurrences(&input);
        let exp: [(&u32, &u32); 5] = [(&3, &4), (&4, &4), (&2, &3), (&1, &2), (&6, &2)];
        let test = exp.iter().zip(obs).all(|(e, o)| *e == o);

        assert_eq!(test, true);
    }

    #[test]
    fn test_convert_to_ranks() {
        let input = vec![(&3, &4), (&4, &4), (&2, &3), (&1, &2), (&6, &2)];

        let obs = _convert_to_ranks(&input);
        let mut exp = HashMap::new();
        exp.insert(4, 4);
        exp.insert(3, 5);
        exp.insert(6, 1);
        exp.insert(2, 3);
        exp.insert(1, 2);
        exp.insert(0, 0);

        assert_eq!(obs.len(), exp.len());
        assert_eq!(exp.keys().all(|k| obs.contains_key(k)), true);
        assert_eq!(exp.iter().all(|(k, v)| obs[k] == *v), true);
    }
}
